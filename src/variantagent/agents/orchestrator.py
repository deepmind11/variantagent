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

    TODO: Implement actual API calls to ClinVar, gnomAD, Ensembl VEP.
    Currently returns an empty annotation as a placeholder.
    """
    start = time.time()
    variant = state["variant"]

    # Placeholder — will be replaced with real API calls
    annotation = VariantAnnotation()
    errors: list[str] = []

    # TODO: Call ClinVar MCP server
    # TODO: Call gnomAD API
    # TODO: Call Ensembl VEP API
    # TODO: Handle rate limiting, retries, fallbacks

    duration_ms = int((time.time() - start) * 1000)
    provenance_entry = ProvenanceEntry(
        step=3,
        agent="annotation_agent",
        action="Database annotation",
        input_summary=f"Variant: {variant.variant_id}",
        output_summary="Queried ClinVar, gnomAD, Ensembl VEP (placeholder)",
        data_source="ClinVar, gnomAD, Ensembl VEP",
        duration_ms=duration_ms,
    )

    return {
        "annotation": annotation,
        "provenance": [provenance_entry],
        "errors": errors,
    }


def literature_node(state: AnalysisState) -> dict[str, Any]:
    """Search PubMed and RAG knowledge base for variant evidence.

    TODO: Implement PubMed search and ACMG guidelines RAG.
    Currently a placeholder.
    """
    start = time.time()
    variant = state["variant"]

    # TODO: PubMed search via NCBI E-utilities
    # TODO: RAG over embedded ACMG guidelines (ChromaDB)

    duration_ms = int((time.time() - start) * 1000)
    provenance_entry = ProvenanceEntry(
        step=4,
        agent="literature_agent",
        action="Literature search",
        input_summary=f"Gene: {variant.gene or 'unknown'}, Variant: {variant.variant_id}",
        output_summary="Literature search (placeholder)",
        data_source="PubMed",
        duration_ms=duration_ms,
    )

    return {
        "provenance": [provenance_entry],
    }


def classification_node(state: AnalysisState) -> dict[str, Any]:
    """Apply ACMG criteria and classify the variant.

    TODO: Implement LLM-based criterion assessment + deterministic rule engine.
    Currently returns a placeholder VUS classification.
    """
    start = time.time()
    variant = state["variant"]

    from variantagent.models.classification import (
        ACMGClassification,
        ACMGClassificationResult,
        ACMGCriteria,
    )

    # Placeholder classification — will be replaced with real ACMG assessment
    classification = ACMGClassification(
        classification=ACMGClassificationResult.VUS,
        criteria=ACMGCriteria(),
        confidence=0.5,
        reasoning="Placeholder — ACMG criterion assessment not yet implemented",
        applied_codes_summary=[],
        classification_rule="No evidence criteria met",
    )

    duration_ms = int((time.time() - start) * 1000)
    provenance_entry = ProvenanceEntry(
        step=5,
        agent="classification_agent",
        action="ACMG classification",
        input_summary=f"Variant: {variant.variant_id}, annotation available: {state['annotation'] is not None}",
        output_summary=f"Classification: {classification.classification.value} "
                       f"(confidence: {classification.confidence})",
        duration_ms=duration_ms,
    )

    return {
        "classification": classification,
        "overall_confidence": classification.confidence,
        "provenance": [provenance_entry],
    }


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
