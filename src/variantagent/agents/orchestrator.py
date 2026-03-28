"""Orchestrator Agent: LangGraph state machine that plans and routes analysis.

This is the central coordination layer. It receives a variant, creates an
analysis plan, routes to specialist agents, handles conditional branching
(QC fail → skip annotation, novel variant → trigger literature), and
assembles the final report with provenance.

The graph structure:

    START → plan → qc → [route] → annotate → [route] → literature → classify
                    ↓ (QC fail)                ↓ (skip lit)            ↓
                  report_warning             classify              review
                                                                     ↓
                                                              [confidence gate]
                                                                ↓          ↓
                                                          hitl_review    report
                                                                ↓
                                                              report → END
"""

from __future__ import annotations

import logging
import operator
import time
import uuid
from typing import Annotated, Any

from langgraph.checkpoint.memory import MemorySaver
from langgraph.graph import END, START, StateGraph
from langgraph.types import Command, interrupt
from typing_extensions import TypedDict

from variantagent.agents.qc_agent import run_qc_assessment
from variantagent.config import settings
from variantagent.models.annotation import VariantAnnotation
from variantagent.models.classification import ACMGClassification
from variantagent.models.qc_metrics import QCAssessment, QCStatus
from variantagent.models.report import ProvenanceEntry, ReviewerFinding, TriageReport
from variantagent.models.variant import Variant

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# State
# ---------------------------------------------------------------------------

class AnalysisState(TypedDict):
    """LangGraph state for the variant analysis workflow.

    Fields without reducers are last-write-wins.
    """

    trace_id: str
    variant: Variant
    sample_id: str | None
    batch_id: str | None

    # Orchestrator planning
    plan: list[str]

    # Agent outputs
    qc_assessment: QCAssessment | None
    annotation: VariantAnnotation | None
    classification: ACMGClassification | None

    # Reviewer
    reviewer_findings: Annotated[list[ReviewerFinding], operator.add]
    overall_confidence: float

    # Report
    report: TriageReport | None

    # Provenance (accumulated across nodes — reducer appends, never overwrites)
    provenance: Annotated[list[ProvenanceEntry], operator.add]
    errors: Annotated[list[str], operator.add]

    # Human-in-the-loop
    requires_human_review: bool
    human_review_reason: str | None


def create_initial_state(
    variant: Variant,
    sample_id: str | None = None,
    batch_id: str | None = None,
) -> AnalysisState:
    """Create initial state for a new analysis run."""
    return AnalysisState(
        trace_id=str(uuid.uuid4()),
        variant=variant,
        sample_id=sample_id,
        batch_id=batch_id,
        plan=[],
        qc_assessment=None,
        annotation=None,
        classification=None,
        reviewer_findings=[],
        overall_confidence=0.0,
        report=None,
        provenance=[],
        errors=[],
        requires_human_review=False,
        human_review_reason=None,
    )


# ---------------------------------------------------------------------------
# Node functions — each receives state, returns partial state update
# ---------------------------------------------------------------------------

def plan_node(state: AnalysisState) -> dict[str, Any]:
    """Create an analysis plan based on the variant and available data."""
    start = time.time()
    variant = state["variant"]

    plan = [
        f"1. QC assessment for sample {state['sample_id'] or 'unknown'}",
        f"2. Query ClinVar, gnomAD, Ensembl VEP for {variant.variant_id}",
    ]

    # Plan literature search if gene is known
    if variant.gene:
        plan.append(f"3. Search PubMed for {variant.gene} variant evidence")
        plan.append(f"4. Apply ACMG criteria for {variant.gene}:{variant.hgvs_p or variant.variant_id}")
    else:
        plan.append("3. Apply ACMG criteria (gene unknown — limited criteria available)")

    plan.append("5. Self-evaluation and contradiction check")
    plan.append("6. Generate structured report with provenance")

    duration_ms = int((time.time() - start) * 1000)
    provenance_entry = ProvenanceEntry(
        step=1,
        agent="orchestrator",
        action="Created analysis plan",
        input_summary=f"Variant: {variant.variant_id}",
        output_summary=f"Plan with {len(plan)} steps",
        duration_ms=duration_ms,
    )

    logger.info("Plan created for %s: %d steps", variant.variant_id, len(plan))

    return {
        "plan": plan,
        "provenance": [provenance_entry],
    }


def qc_node(state: AnalysisState) -> dict[str, Any]:
    """Run QC assessment on the sample."""
    start = time.time()
    sample_id = state["sample_id"] or "unknown"

    # In a real system, we'd load QC files here (flagstat, MultiQC).
    # For now, we run the assessment with whatever data is available.
    # The QC agent handles the case where no data is present (passes by default).
    qc_assessment = run_qc_assessment(
        sample_id=sample_id,
        flagstat=None,  # TODO: load from file path when available
        multiqc=None,   # TODO: load from file path when available
        variant_region_coverage=state["variant"].depth,
    )

    duration_ms = int((time.time() - start) * 1000)
    provenance_entry = ProvenanceEntry(
        step=2,
        agent="qc_agent",
        action="QC assessment",
        input_summary=f"Sample: {sample_id}",
        output_summary=f"Status: {qc_assessment.overall_status.value}, "
                       f"{len(qc_assessment.issues)} issues",
        duration_ms=duration_ms,
    )

    logger.info(
        "QC for %s: %s (%d issues)",
        sample_id,
        qc_assessment.overall_status.value,
        len(qc_assessment.issues),
    )

    return {
        "qc_assessment": qc_assessment,
        "provenance": [provenance_entry],
    }


def annotation_node(state: AnalysisState) -> dict[str, Any]:
    """Query public databases for variant annotation.

    Runs ClinVar + Ensembl VEP in parallel (both are fast),
    then gnomAD sequentially (aggressive rate limiting).
    If any API fails, continues with the others (graceful degradation).
    """
    import asyncio

    start = time.time()
    variant = state["variant"]

    annotation, errors, sources = _run_annotation_sync(variant)

    duration_ms = int((time.time() - start) * 1000)

    # Build summary of what was found
    found_sources = []
    if annotation.clinvar.found:
        found_sources.append(f"ClinVar: {annotation.clinvar.clinical_significance}")
    if annotation.ensembl_vep.found:
        found_sources.append(f"VEP: {annotation.ensembl_vep.consequence_type}")
    if annotation.gnomad.found:
        af_str = f"{annotation.gnomad.overall_af:.6f}" if annotation.gnomad.overall_af else "N/A"
        found_sources.append(f"gnomAD: AF={af_str}")

    output_summary = "; ".join(found_sources) if found_sources else "No databases returned results"

    provenance_entry = ProvenanceEntry(
        step=3,
        agent="annotation_agent",
        action="Database annotation",
        input_summary=f"Variant: {variant.variant_id}",
        output_summary=output_summary,
        data_source=", ".join(sources),
        duration_ms=duration_ms,
    )

    logger.info("Annotation for %s: %s", variant.variant_id, output_summary)

    return {
        "annotation": annotation,
        "provenance": [provenance_entry],
        "errors": errors,
    }


def _run_annotation_sync(variant: Variant) -> tuple[VariantAnnotation, list[str], list[str]]:
    """Run all annotation queries synchronously.

    Uses synchronous httpx to avoid event loop conflicts with LangGraph.
    ClinVar + VEP run sequentially (fast enough), gnomAD last (slow due to rate limiting).

    Returns:
        Tuple of (annotation, errors, sources_queried).
    """
    import httpx

    from variantagent.tools.clinvar_client import ClinVarAnnotation, _build_query, _esearch, _esummary, _parse_esummary
    from variantagent.models.annotation import EnsemblVEPAnnotation, GnomADFrequency

    errors: list[str] = []
    sources: list[str] = ["ClinVar", "Ensembl VEP", "gnomAD"]

    clinvar_result = ClinVarAnnotation(found=False)
    vep_result = EnsemblVEPAnnotation(found=False)
    gnomad_result = GnomADFrequency(found=False)

    with httpx.Client(timeout=30.0) as client:
        # ClinVar (sync)
        try:
            from variantagent.tools.clinvar_client import _base_params, EUTILS_BASE

            query = _build_query(variant)
            search_params = {
                **_base_params(),
                "db": "clinvar", "term": query, "retmode": "json", "retmax": "5",
            }
            resp = client.get(f"{EUTILS_BASE}/esearch.fcgi", params=search_params)
            resp.raise_for_status()
            uids = resp.json().get("esearchresult", {}).get("idlist", [])

            if uids:
                summary_params = {
                    **_base_params(),
                    "db": "clinvar", "id": ",".join(uids), "retmode": "json",
                }
                resp = client.get(f"{EUTILS_BASE}/esummary.fcgi", params=summary_params)
                resp.raise_for_status()
                clinvar_result = _parse_esummary(resp.json(), uids)
        except Exception as e:
            errors.append(f"ClinVar error: {e}")
            logger.error("ClinVar sync query failed: %s", e)

        # Ensembl VEP (sync)
        try:
            from variantagent.tools.ensembl_client import _build_vep_url, _parse_vep_response

            url = _build_vep_url(variant)
            params = {
                "content-type": "application/json",
                "pick": "1", "SIFT": "b", "PolyPhen": "b",
                "domains": "1", "canonical": "1", "hgvs": "1",
            }
            resp = client.get(url, params=params, headers={"Content-Type": "application/json"})
            if resp.status_code != 429:
                resp.raise_for_status()
                vep_result = _parse_vep_response(resp.json())
        except Exception as e:
            errors.append(f"VEP error: {e}")
            logger.error("VEP sync query failed: %s", e)

        # gnomAD (sync — GraphQL POST)
        try:
            from variantagent.tools.gnomad_client import (
                GNOMAD_API, VARIANT_QUERY, _build_variant_id, _parse_gnomad_response,
            )

            payload = {
                "query": VARIANT_QUERY,
                "variables": {"variantId": _build_variant_id(variant), "datasetId": "gnomad_r4"},
            }
            resp = client.post(
                GNOMAD_API, json=payload,
                headers={"Content-Type": "application/json"},
            )
            resp.raise_for_status()
            data = resp.json()
            if "errors" not in data:
                gnomad_result = _parse_gnomad_response(data)
        except Exception as e:
            errors.append(f"gnomAD error: {e}")
            logger.error("gnomAD sync query failed: %s", e)

    annotation = VariantAnnotation(
        clinvar=clinvar_result,
        gnomad=gnomad_result,
        ensembl_vep=vep_result,
        annotation_errors=errors,
    )

    return annotation, errors, sources


def literature_node(state: AnalysisState) -> dict[str, Any]:
    """Search PubMed for variant-specific literature evidence.

    Uses progressive broadening: gene + HGVS → gene + protein change → gene + disease.
    RAG over ACMG guidelines will be added in a future step.
    """
    import asyncio

    start = time.time()
    variant = state["variant"]

    pmids: list[str] = []
    errors: list[str] = []

    if variant.gene:
        articles, pub_err = _run_literature_search_sync(variant)
        if pub_err:
            errors.append(pub_err)
        pmids = [a.pmid for a in articles]

    duration_ms = int((time.time() - start) * 1000)

    output_summary = f"Found {len(pmids)} articles" if pmids else "No articles found"
    if not variant.gene:
        output_summary = "Skipped — no gene symbol available"

    provenance_entry = ProvenanceEntry(
        step=4,
        agent="literature_agent",
        action="Literature search",
        input_summary=f"Gene: {variant.gene or 'unknown'}, Variant: {variant.variant_id}",
        output_summary=output_summary,
        data_source="PubMed",
        duration_ms=duration_ms,
    )

    # Store PMIDs in the annotation if available
    result: dict[str, Any] = {
        "provenance": [provenance_entry],
        "errors": errors,
    }

    # Update annotation with PubMed references
    if pmids and state["annotation"]:
        updated_annotation = state["annotation"].model_copy(
            update={"pubmed_references": pmids}
        )
        result["annotation"] = updated_annotation

    return result


def _run_literature_search_sync(variant: Variant) -> tuple[list, str | None]:
    """Run PubMed search synchronously."""
    import httpx

    from variantagent.tools.pubmed_client import PubMedArticle, _build_search_queries
    from variantagent.tools.clinvar_client import _base_params, EUTILS_BASE

    if not variant.gene:
        return [], None

    queries = _build_search_queries(variant.gene, variant.hgvs_p, variant.variant_id)
    articles: list[PubMedArticle] = []
    seen: set[str] = set()

    try:
        with httpx.Client(timeout=30.0) as client:
            all_pmids: list[str] = []

            for query in queries:
                if len(all_pmids) >= 5:
                    break
                params = {
                    **_base_params(),
                    "db": "pubmed", "term": query, "retmode": "json",
                    "retmax": "5", "sort": "relevance",
                }
                resp = client.get(f"{EUTILS_BASE}/esearch.fcgi", params=params)
                resp.raise_for_status()
                pmids = resp.json().get("esearchresult", {}).get("idlist", [])
                for pmid in pmids:
                    if pmid not in seen and len(all_pmids) < 5:
                        all_pmids.append(pmid)
                        seen.add(pmid)

            if all_pmids:
                params = {
                    **_base_params(),
                    "db": "pubmed", "id": ",".join(all_pmids), "retmode": "json",
                }
                resp = client.get(f"{EUTILS_BASE}/esummary.fcgi", params=params)
                resp.raise_for_status()
                result = resp.json().get("result", {})

                for pmid in all_pmids:
                    record = result.get(pmid, {})
                    if record and "error" not in record:
                        authors = [a.get("name", "") for a in record.get("authors", [])]
                        articles.append(PubMedArticle(
                            pmid=pmid,
                            title=record.get("title", ""),
                            journal=record.get("fulljournalname", ""),
                            year=record.get("pubdate", "")[:4],
                            authors=authors,
                        ))

        return articles, None
    except Exception as e:
        logger.error("PubMed sync search failed: %s", e)
        return [], f"PubMed error: {e}"


def classification_node(state: AnalysisState) -> dict[str, Any]:
    """Apply ACMG criteria and classify the variant.

    Two-phase approach:
    1. Rule-based criterion assessment using annotation data (no LLM needed
       for criteria that can be evaluated deterministically from database evidence)
    2. Deterministic ACMG combining rules via acmg_engine.classify()

    LLM-based criterion assessment for ambiguous criteria will be added
    once the rule-based approach is validated.
    """
    start = time.time()
    variant = state["variant"]
    annotation = state["annotation"]

    from variantagent.models.classification import (
        ACMGClassification,
        ACMGCriteria,
        EvidenceCode,
        EvidenceDirection,
        EvidenceStrength,
    )
    from variantagent.tools.acmg_engine import classify as acmg_classify

    criteria = ACMGCriteria()

    # Evaluate criteria from annotation evidence (deterministic, no LLM)
    if annotation:
        criteria = _evaluate_criteria_from_evidence(variant, annotation)

    # Run deterministic ACMG combining rules
    classification_result, rule_description = acmg_classify(criteria)

    # Calculate confidence based on evidence completeness
    confidence = _calculate_confidence(annotation, criteria)

    applied_codes = criteria.get_applied_codes()
    applied_summary = [c.code for c in applied_codes]

    reasoning_parts = [f"Classification: {classification_result.value}"]
    reasoning_parts.append(f"Rule applied: {rule_description}")
    if applied_codes:
        reasoning_parts.append(f"Evidence codes: {', '.join(applied_summary)}")
        for code in applied_codes:
            reasoning_parts.append(f"  {code.code}: {code.reasoning}")
    else:
        reasoning_parts.append("No ACMG evidence criteria were met from available data.")

    classification = ACMGClassification(
        classification=classification_result,
        criteria=criteria,
        confidence=confidence,
        reasoning="\n".join(reasoning_parts),
        applied_codes_summary=applied_summary,
        classification_rule=rule_description,
    )

    duration_ms = int((time.time() - start) * 1000)
    provenance_entry = ProvenanceEntry(
        step=5,
        agent="classification_agent",
        action="ACMG classification",
        input_summary=f"Variant: {variant.variant_id}, "
                      f"ClinVar: {annotation.clinvar.found if annotation else 'N/A'}, "
                      f"gnomAD AF: {annotation.gnomad.overall_af if annotation and annotation.gnomad.found else 'N/A'}",
        output_summary=f"{classification_result.value} ({', '.join(applied_summary) or 'no criteria met'}) "
                       f"confidence: {confidence:.2f}",
        duration_ms=duration_ms,
    )

    logger.info(
        "Classification for %s: %s (confidence: %.2f, codes: %s)",
        variant.variant_id,
        classification_result.value,
        confidence,
        applied_summary,
    )

    return {
        "classification": classification,
        "overall_confidence": confidence,
        "provenance": [provenance_entry],
    }


def _evaluate_criteria_from_evidence(
    variant: Variant,
    annotation: VariantAnnotation,
) -> "ACMGCriteria":
    """Evaluate ACMG criteria deterministically from annotation evidence.

    These criteria can be assessed without LLM judgment — they depend on
    objective data from databases.
    """
    from variantagent.models.classification import (
        ACMGCriteria,
        EvidenceCode,
        EvidenceDirection,
        EvidenceStrength,
    )

    kwargs: dict[str, EvidenceCode] = {}

    # --- BA1: Allele frequency > 5% (standalone benign) ---
    if annotation.gnomad.found and annotation.gnomad.overall_af is not None:
        if annotation.gnomad.overall_af > 0.05:
            kwargs["ba1"] = EvidenceCode(
                code="BA1", name="Allele frequency > 5%",
                direction=EvidenceDirection.BENIGN, strength=EvidenceStrength.VERY_STRONG,
                applied=True,
                reasoning=f"gnomAD overall AF = {annotation.gnomad.overall_af:.4f} (> 0.05 threshold)",
                data_source="gnomAD", confidence=0.99,
            )

    # --- BS1: Allele frequency greater than expected for disorder ---
    if annotation.gnomad.found and annotation.gnomad.overall_af is not None:
        if 0.01 < annotation.gnomad.overall_af <= 0.05:
            kwargs["bs1"] = EvidenceCode(
                code="BS1", name="Allele frequency greater than expected",
                direction=EvidenceDirection.BENIGN, strength=EvidenceStrength.STRONG,
                applied=True,
                reasoning=f"gnomAD overall AF = {annotation.gnomad.overall_af:.4f} (> 0.01, suggesting common variant)",
                data_source="gnomAD", confidence=0.90,
            )

    # --- PM2: Absent from population databases ---
    if annotation.gnomad.found and annotation.gnomad.overall_af is not None:
        if annotation.gnomad.overall_af < 0.0001:
            kwargs["pm2"] = EvidenceCode(
                code="PM2", name="Absent/extremely low frequency in population databases",
                direction=EvidenceDirection.PATHOGENIC, strength=EvidenceStrength.MODERATE,
                applied=True,
                reasoning=f"gnomAD overall AF = {annotation.gnomad.overall_af:.6f} (< 0.0001 threshold)",
                data_source="gnomAD", confidence=0.85,
            )
    elif annotation.gnomad.found is False:
        kwargs["pm2"] = EvidenceCode(
            code="PM2", name="Absent from population databases",
            direction=EvidenceDirection.PATHOGENIC, strength=EvidenceStrength.MODERATE,
            applied=True,
            reasoning="Variant not found in gnomAD — absent from population databases",
            data_source="gnomAD", confidence=0.80,
        )

    # --- PP5: Reputable source reports pathogenic ---
    if annotation.clinvar.found and annotation.clinvar.clinical_significance:
        sig = annotation.clinvar.clinical_significance.lower()
        stars = annotation.clinvar.review_stars or 0

        if "pathogenic" in sig and "likely" not in sig and stars >= 2:
            kwargs["pp5"] = EvidenceCode(
                code="PP5", name="Reputable source reports pathogenic",
                direction=EvidenceDirection.PATHOGENIC, strength=EvidenceStrength.SUPPORTING,
                applied=True,
                reasoning=f"ClinVar: '{annotation.clinvar.clinical_significance}' "
                          f"({stars}-star review, {annotation.clinvar.submitter_count} submitters)",
                data_source="ClinVar", confidence=0.90,
            )
        elif "likely pathogenic" in sig and stars >= 2:
            kwargs["pp5"] = EvidenceCode(
                code="PP5", name="Reputable source reports likely pathogenic",
                direction=EvidenceDirection.PATHOGENIC, strength=EvidenceStrength.SUPPORTING,
                applied=True,
                reasoning=f"ClinVar: '{annotation.clinvar.clinical_significance}' ({stars}-star review)",
                data_source="ClinVar", confidence=0.80,
            )

    # --- BP6: Reputable source reports benign ---
    if annotation.clinvar.found and annotation.clinvar.clinical_significance:
        sig = annotation.clinvar.clinical_significance.lower()
        stars = annotation.clinvar.review_stars or 0

        if ("benign" in sig or "likely benign" in sig) and stars >= 2:
            kwargs["bp6"] = EvidenceCode(
                code="BP6", name="Reputable source reports benign",
                direction=EvidenceDirection.BENIGN, strength=EvidenceStrength.SUPPORTING,
                applied=True,
                reasoning=f"ClinVar: '{annotation.clinvar.clinical_significance}' ({stars}-star review)",
                data_source="ClinVar", confidence=0.90,
            )

    # --- PP3: Computational evidence supports deleterious ---
    if annotation.ensembl_vep.found:
        sift_del = annotation.ensembl_vep.sift_prediction == "deleterious"
        polyphen_dam = annotation.ensembl_vep.polyphen_prediction in (
            "probably_damaging", "possibly_damaging"
        )
        if sift_del and polyphen_dam:
            kwargs["pp3"] = EvidenceCode(
                code="PP3", name="Computational evidence supports deleterious",
                direction=EvidenceDirection.PATHOGENIC, strength=EvidenceStrength.SUPPORTING,
                applied=True,
                reasoning=f"SIFT: {annotation.ensembl_vep.sift_prediction} "
                          f"(score: {annotation.ensembl_vep.sift_score}), "
                          f"PolyPhen: {annotation.ensembl_vep.polyphen_prediction} "
                          f"(score: {annotation.ensembl_vep.polyphen_score})",
                data_source="Ensembl VEP", confidence=0.70,
            )

    # --- BP4: Computational evidence supports benign ---
    if annotation.ensembl_vep.found:
        sift_tol = annotation.ensembl_vep.sift_prediction == "tolerated"
        polyphen_ben = annotation.ensembl_vep.polyphen_prediction == "benign"
        if sift_tol and polyphen_ben:
            kwargs["bp4"] = EvidenceCode(
                code="BP4", name="Computational evidence supports benign",
                direction=EvidenceDirection.BENIGN, strength=EvidenceStrength.SUPPORTING,
                applied=True,
                reasoning=f"SIFT: tolerated (score: {annotation.ensembl_vep.sift_score}), "
                          f"PolyPhen: benign (score: {annotation.ensembl_vep.polyphen_score})",
                data_source="Ensembl VEP", confidence=0.70,
            )

    # --- PM1: Located in mutational hot spot / functional domain ---
    if annotation.ensembl_vep.found and annotation.ensembl_vep.protein_domain:
        kwargs["pm1"] = EvidenceCode(
            code="PM1", name="Located in mutational hot spot / functional domain",
            direction=EvidenceDirection.PATHOGENIC, strength=EvidenceStrength.MODERATE,
            applied=True,
            reasoning=f"Variant falls in protein domain: {annotation.ensembl_vep.protein_domain}",
            data_source="Ensembl VEP", confidence=0.65,
        )

    return ACMGCriteria(**kwargs)


def _calculate_confidence(
    annotation: VariantAnnotation | None,
    criteria: "ACMGCriteria",
) -> float:
    """Calculate confidence score based on evidence completeness.

    Higher confidence when more data sources contributed evidence.
    """
    if annotation is None:
        return 0.2

    score = 0.3  # Base confidence

    # Bonus for each data source that returned results
    if annotation.clinvar.found:
        score += 0.15
    if annotation.gnomad.found:
        score += 0.15
    if annotation.ensembl_vep.found:
        score += 0.10

    # Bonus for number of criteria evaluated
    applied = criteria.get_applied_codes()
    score += min(len(applied) * 0.05, 0.20)

    # Bonus for ClinVar review quality
    if annotation.clinvar.found and (annotation.clinvar.review_stars or 0) >= 2:
        score += 0.10

    return min(score, 1.0)


def review_node(state: AnalysisState) -> dict[str, Any]:
    """Self-evaluation: cross-check conclusions and detect contradictions.

    TODO: Implement claim extraction, source verification, contradiction detection.
    Currently performs basic consistency checks.
    """
    start = time.time()
    findings: list[ReviewerFinding] = []

    # Basic consistency check: QC said unreliable but we still classified
    qc = state["qc_assessment"]
    classification = state["classification"]

    if qc and not qc.reliable_for_interpretation and classification:
        findings.append(
            ReviewerFinding(
                claim=f"Variant classified as {classification.classification.value}",
                supported=False,
                concern="QC assessment indicates variant call may not be reliable, "
                        "but classification was still performed. Interpret with extreme caution.",
                hallucination_risk="high",
            )
        )

    # Check if classification has low confidence
    if classification and classification.confidence < 0.5:
        findings.append(
            ReviewerFinding(
                claim=f"Classification confidence: {classification.confidence}",
                supported=True,
                concern="Low confidence classification — insufficient evidence to "
                        "support a definitive call.",
                hallucination_risk="medium",
            )
        )

    # Determine if human review is needed
    confidence = state["overall_confidence"]
    threshold = settings.hitl_confidence_threshold
    needs_review = confidence < threshold

    review_reason = None
    if needs_review:
        review_reason = (
            f"Overall confidence ({confidence:.2f}) is below threshold ({threshold:.2f}). "
            f"{len(findings)} reviewer concern(s) raised."
        )

    duration_ms = int((time.time() - start) * 1000)
    provenance_entry = ProvenanceEntry(
        step=6,
        agent="reviewer_agent",
        action="Self-evaluation",
        input_summary="Cross-checking QC, annotation, and classification consistency",
        output_summary=f"{len(findings)} findings, needs_review={needs_review}",
        duration_ms=duration_ms,
    )

    return {
        "reviewer_findings": findings,
        "requires_human_review": needs_review,
        "human_review_reason": review_reason,
        "provenance": [provenance_entry],
    }


def hitl_node(state: AnalysisState) -> dict[str, Any]:
    """Human-in-the-loop checkpoint.

    Pauses execution and presents findings for human review.
    Resumes when the human provides a decision via Command(resume=...).
    """
    human_decision = interrupt({
        "variant": state["variant"].variant_id,
        "classification": state["classification"].classification.value if state["classification"] else "unknown",
        "confidence": state["overall_confidence"],
        "reason": state["human_review_reason"],
        "reviewer_concerns": [f.concern for f in state["reviewer_findings"] if f.concern],
        "prompt": "Review the classification and concerns above. "
                  "Respond with {'approve': true/false, 'override_classification': '...' (optional)}",
    })

    # Process human decision
    approved = human_decision.get("approve", False) if isinstance(human_decision, dict) else False
    override = human_decision.get("override_classification") if isinstance(human_decision, dict) else None

    provenance_entry = ProvenanceEntry(
        step=7,
        agent="human",
        action="Human review checkpoint",
        input_summary=f"Confidence: {state['overall_confidence']:.2f}, "
                      f"{len(state['reviewer_findings'])} concerns",
        output_summary=f"Approved: {approved}, Override: {override or 'none'}",
    )

    result: dict[str, Any] = {
        "requires_human_review": False,
        "provenance": [provenance_entry],
    }

    # If human provided an override classification, update it
    if override and state["classification"]:
        from variantagent.models.classification import ACMGClassificationResult

        try:
            new_class = ACMGClassificationResult(override)
            updated = state["classification"].model_copy(
                update={"classification": new_class, "reasoning": f"Human override: {override}"}
            )
            result["classification"] = updated
        except ValueError:
            result["errors"] = [f"Invalid override classification: {override}"]

    return result


def report_node(state: AnalysisState) -> dict[str, Any]:
    """Generate the final TriageReport with full provenance."""
    start = time.time()

    report = TriageReport(
        trace_id=state["trace_id"],
        variant=state["variant"],
        sample_id=state["sample_id"],
        batch_id=state["batch_id"],
        qc_assessment=state["qc_assessment"],
        annotation=state["annotation"],
        classification=state["classification"],
        reviewer_findings=state["reviewer_findings"],
        overall_confidence=state["overall_confidence"],
        requires_human_review=state["requires_human_review"],
        human_review_reason=state["human_review_reason"],
        provenance=state["provenance"],
        analysis_plan=state["plan"],
        limitations=[
            "ACMG criterion assessment not yet implemented — placeholder classification",
            "Database annotations are placeholders — no live API calls yet",
            "Literature search not yet implemented",
            "Batch comparison not yet implemented",
            "This is a development scaffold, not for clinical use",
        ],
    )

    # Generate natural language summary
    qc_text = "not assessed"
    if state["qc_assessment"]:
        qc_text = state["qc_assessment"].overall_status.value

    qc_aborted = (
        state["qc_assessment"] is not None
        and not state["qc_assessment"].reliable_for_interpretation
        and state["classification"] is None
    )

    if qc_aborted:
        report.natural_language_summary = (
            f"Variant {state['variant'].variant_id} "
            f"(gene: {state['variant'].gene or 'unknown'}) "
            f"was NOT assessed due to QC failure (status: {qc_text}). "
            f"Variant interpretation was skipped because the sequencing data "
            f"is unreliable. Recommended action: {state['qc_assessment'].issues[0].recommended_action if state['qc_assessment'].issues else 'review QC metrics'}."
        )
    else:
        classification_text = "Uncertain Significance"
        if state["classification"]:
            classification_text = state["classification"].classification.value

        report.natural_language_summary = (
            f"Variant {state['variant'].variant_id} "
            f"(gene: {state['variant'].gene or 'unknown'}) "
            f"was classified as {classification_text} "
            f"with {state['overall_confidence']:.0%} confidence. "
            f"QC status: {qc_text}. "
            f"{len(state['reviewer_findings'])} reviewer concern(s)."
        )

    duration_ms = int((time.time() - start) * 1000)

    logger.info(
        "Report generated for %s: %s (confidence: %.2f)",
        state["variant"].variant_id,
        classification_text,
        state["overall_confidence"],
    )

    return {"report": report}


# ---------------------------------------------------------------------------
# Routing functions — decide which node to go to next
# ---------------------------------------------------------------------------

def route_after_qc(state: AnalysisState) -> str:
    """Decide what to do after QC assessment.

    If QC failed and the variant call is unreliable, skip straight to
    the report (with a warning). Otherwise, continue to annotation.
    """
    qc = state["qc_assessment"]
    if qc and qc.overall_status == QCStatus.FAIL and not qc.reliable_for_interpretation:
        logger.warning("QC failed — skipping annotation, going to report")
        return "report"
    return "annotate"


def route_after_annotation(state: AnalysisState) -> str:
    """Decide whether to search literature before classification.

    If the variant is novel (not found in ClinVar) and a gene is known,
    trigger literature search. Otherwise, skip to classification.
    """
    annotation = state["annotation"]
    variant = state["variant"]

    # Trigger literature search if: variant has a known gene AND
    # ClinVar didn't find it (novel variant needs more evidence)
    if variant.gene and annotation and not annotation.clinvar.found:
        logger.info("Novel variant in %s — triggering literature search", variant.gene)
        return "literature"

    return "classify"


def route_after_review(state: AnalysisState) -> str:
    """Decide whether human review is needed after self-evaluation."""
    if state["requires_human_review"]:
        return "hitl"
    return "report"


# ---------------------------------------------------------------------------
# Graph construction
# ---------------------------------------------------------------------------

def build_graph(checkpointer: Any | None = None) -> StateGraph:
    """Build and compile the VariantAgent LangGraph workflow.

    Args:
        checkpointer: LangGraph checkpointer for state persistence.
            Required for human-in-the-loop. Defaults to MemorySaver
            if None.

    Returns:
        Compiled LangGraph StateGraph ready for invocation.
    """
    if checkpointer is None:
        checkpointer = MemorySaver()

    builder = StateGraph(AnalysisState)

    # Add nodes
    builder.add_node("plan", plan_node)
    builder.add_node("qc", qc_node)
    builder.add_node("annotate", annotation_node)
    builder.add_node("literature", literature_node)
    builder.add_node("classify", classification_node)
    builder.add_node("review", review_node)
    builder.add_node("hitl", hitl_node)
    builder.add_node("report", report_node)

    # Edges
    builder.add_edge(START, "plan")
    builder.add_edge("plan", "qc")

    # After QC: continue or skip to report
    builder.add_conditional_edges(
        "qc",
        route_after_qc,
        {"annotate": "annotate", "report": "report"},
    )

    # After annotation: literature search or straight to classification
    builder.add_conditional_edges(
        "annotate",
        route_after_annotation,
        {"literature": "literature", "classify": "classify"},
    )

    builder.add_edge("literature", "classify")
    builder.add_edge("classify", "review")

    # After review: human checkpoint or report
    builder.add_conditional_edges(
        "review",
        route_after_review,
        {"hitl": "hitl", "report": "report"},
    )

    builder.add_edge("hitl", "report")
    builder.add_edge("report", END)

    return builder.compile(checkpointer=checkpointer)


# ---------------------------------------------------------------------------
# Convenience runner
# ---------------------------------------------------------------------------

def analyze_variant(
    variant: Variant,
    sample_id: str | None = None,
    batch_id: str | None = None,
    thread_id: str | None = None,
    auto_approve: bool = True,
) -> TriageReport:
    """Run the full variant analysis pipeline.

    This is the main entry point for programmatic use.

    Args:
        variant: The variant to analyze.
        sample_id: Optional sample identifier.
        batch_id: Optional batch identifier.
        thread_id: Optional thread ID for checkpointing (auto-generated if None).
        auto_approve: If True, automatically approve at HITL checkpoints
            instead of blocking. Set to False for interactive use where
            a human will provide input via Command(resume=...).

    Returns:
        Complete TriageReport with provenance.

    Raises:
        RuntimeError: If the graph fails to produce a report.
    """
    graph = build_graph()
    initial_state = create_initial_state(variant, sample_id, batch_id)

    config = {"configurable": {"thread_id": thread_id or initial_state["trace_id"]}}

    result = graph.invoke(initial_state, config=config)

    # If the graph hit a HITL interrupt, auto-approve and resume
    if auto_approve and "__interrupt__" in result:
        logger.info("HITL interrupt hit — auto-approving (auto_approve=True)")
        result = graph.invoke(
            Command(resume={"approve": True}),
            config=config,
        )

    report = result.get("report")
    if not isinstance(report, TriageReport):
        raise RuntimeError(
            f"Graph did not produce a TriageReport. Final state keys: {list(result.keys())}"
        )

    return report
