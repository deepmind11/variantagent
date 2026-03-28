"""Integration tests for the Orchestrator LangGraph workflow.

These tests mock the external API calls to run without network access.
"""

from unittest.mock import MagicMock, patch

import pytest

from variantagent.agents.orchestrator import analyze_variant
from variantagent.models.annotation import ClinVarAnnotation, EnsemblVEPAnnotation, GnomADFrequency, VariantAnnotation
from variantagent.models.qc_metrics import QCStatus
from variantagent.models.report import TriageReport
from variantagent.models.variant import Variant, VariantType


def _mock_annotation_sync(variant):
    """Return a realistic mock annotation — variant not found in ClinVar."""
    annotation = VariantAnnotation(
        clinvar=ClinVarAnnotation(found=False),
        gnomad=GnomADFrequency(found=True, overall_af=0.0001, allele_count=10, allele_number=100000),
        ensembl_vep=EnsemblVEPAnnotation(
            found=True, consequence_type="missense_variant", impact="MODERATE",
            gene_symbol="TP53",
        ),
    )
    return annotation, [], ["ClinVar", "Ensembl VEP", "gnomAD"]


def _mock_annotation_sync_clinvar_found(variant):
    """Return a mock annotation where ClinVar DOES find the variant."""
    annotation = VariantAnnotation(
        clinvar=ClinVarAnnotation(
            found=True, clinical_significance="Pathogenic",
            review_status="criteria provided, multiple submitters, no conflicts",
            review_stars=2,
        ),
        gnomad=GnomADFrequency(found=True, overall_af=0.00001),
        ensembl_vep=EnsemblVEPAnnotation(found=True, consequence_type="missense_variant", impact="MODERATE"),
    )
    return annotation, [], ["ClinVar", "Ensembl VEP", "gnomAD"]


def _mock_literature_sync(variant):
    """Return empty literature results."""
    return [], None


@pytest.fixture(autouse=True)
def mock_external_apis():
    """Mock all external API calls for integration tests."""
    with (
        patch(
            "variantagent.agents.orchestrator._run_annotation_sync",
            side_effect=_mock_annotation_sync,
        ),
        patch(
            "variantagent.agents.orchestrator._run_literature_search_sync",
            side_effect=_mock_literature_sync,
        ),
    ):
        yield


class TestOrchestratorGraph:
    """Test the full LangGraph workflow end-to-end."""

    def test_full_pipeline_produces_report(self, tp53_missense: Variant) -> None:
        report = analyze_variant(tp53_missense, sample_id="S001")
        assert isinstance(report, TriageReport)
        assert report.variant.variant_id == tp53_missense.variant_id
        assert report.sample_id == "S001"
        assert report.trace_id != ""

    def test_report_has_plan(self, tp53_missense: Variant) -> None:
        report = analyze_variant(tp53_missense)
        assert len(report.analysis_plan) > 0
        assert any("QC" in step for step in report.analysis_plan)

    def test_report_has_provenance(self, tp53_missense: Variant) -> None:
        report = analyze_variant(tp53_missense)
        assert len(report.provenance) >= 3
        agents_involved = {p.agent for p in report.provenance}
        assert "orchestrator" in agents_involved
        assert "qc_agent" in agents_involved

    def test_report_has_classification(self, tp53_missense: Variant) -> None:
        report = analyze_variant(tp53_missense)
        assert report.classification is not None

    def test_report_has_limitations(self, tp53_missense: Variant) -> None:
        report = analyze_variant(tp53_missense)
        assert len(report.limitations) > 0

    def test_report_has_natural_language_summary(self, tp53_missense: Variant) -> None:
        report = analyze_variant(tp53_missense)
        assert report.natural_language_summary != ""

    def test_variant_without_gene(self) -> None:
        variant = Variant(
            chromosome="chr1", position=100000, reference="A", alternate="G",
            variant_type=VariantType.SNV,
        )
        report = analyze_variant(variant)
        assert isinstance(report, TriageReport)
        assert report.classification is not None

    def test_different_variants_get_different_trace_ids(
        self, tp53_missense: Variant, brca1_frameshift: Variant
    ) -> None:
        report1 = analyze_variant(tp53_missense)
        report2 = analyze_variant(brca1_frameshift)
        assert report1.trace_id != report2.trace_id


class TestDynamicRouting:
    """Test that conditional edges route correctly."""

    def test_qc_pass_continues_to_annotation(self, tp53_missense: Variant) -> None:
        report = analyze_variant(tp53_missense)
        agent_sequence = [p.agent for p in report.provenance]
        assert "annotation_agent" in agent_sequence

    def test_novel_variant_triggers_literature(self, tp53_missense: Variant) -> None:
        """ClinVar found=False + gene present → literature search triggered."""
        report = analyze_variant(tp53_missense)
        agent_sequence = [p.agent for p in report.provenance]
        assert "literature_agent" in agent_sequence

    def test_clinvar_found_skips_literature(self, tp53_missense: Variant) -> None:
        """ClinVar found=True → literature search skipped, goes straight to classification."""
        with patch(
            "variantagent.agents.orchestrator._run_annotation_sync",
            side_effect=_mock_annotation_sync_clinvar_found,
        ):
            report = analyze_variant(tp53_missense)
            agent_sequence = [p.agent for p in report.provenance]
            # Literature should NOT be in the sequence when ClinVar finds the variant
            assert "literature_agent" not in agent_sequence


class TestAnnotationIntegration:
    """Test that annotation data flows through to the report."""

    def test_annotation_present_in_report(self, tp53_missense: Variant) -> None:
        report = analyze_variant(tp53_missense)
        assert report.annotation is not None
        assert report.annotation.ensembl_vep.found is True
        assert report.annotation.gnomad.found is True

    def test_provenance_shows_annotation_details(self, tp53_missense: Variant) -> None:
        report = analyze_variant(tp53_missense)
        annotation_provenance = [p for p in report.provenance if p.agent == "annotation_agent"]
        assert len(annotation_provenance) == 1
        assert "gnomAD" in annotation_provenance[0].data_source
