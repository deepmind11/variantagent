"""Integration tests for the Orchestrator LangGraph workflow."""

import pytest

from variantagent.agents.orchestrator import (
    AnalysisState,
    analyze_variant,
    build_graph,
    create_initial_state,
)
from variantagent.models.qc_metrics import QCStatus
from variantagent.models.report import TriageReport
from variantagent.models.variant import Variant, VariantType


class TestOrchestratorGraph:
    """Test the full LangGraph workflow end-to-end."""

    def test_full_pipeline_produces_report(self, tp53_missense: Variant) -> None:
        """The graph should produce a TriageReport for a standard variant."""
        report = analyze_variant(tp53_missense, sample_id="S001")

        assert isinstance(report, TriageReport)
        assert report.variant.variant_id == tp53_missense.variant_id
        assert report.sample_id == "S001"
        assert report.trace_id != ""

    def test_report_has_plan(self, tp53_missense: Variant) -> None:
        """The report should contain the analysis plan."""
        report = analyze_variant(tp53_missense)

        assert len(report.analysis_plan) > 0
        assert any("QC" in step for step in report.analysis_plan)

    def test_report_has_provenance(self, tp53_missense: Variant) -> None:
        """Every step should be tracked in provenance."""
        report = analyze_variant(tp53_missense)

        assert len(report.provenance) >= 3  # plan + qc + at least one more
        agents_involved = {p.agent for p in report.provenance}
        assert "orchestrator" in agents_involved
        assert "qc_agent" in agents_involved

    def test_report_has_classification(self, tp53_missense: Variant) -> None:
        """The report should contain a classification (even if placeholder)."""
        report = analyze_variant(tp53_missense)

        assert report.classification is not None
        assert report.classification.classification is not None

    def test_report_has_limitations(self, tp53_missense: Variant) -> None:
        """The report should honestly state its limitations."""
        report = analyze_variant(tp53_missense)

        assert len(report.limitations) > 0
        assert any("placeholder" in l.lower() or "not yet" in l.lower()
                    for l in report.limitations)

    def test_report_has_natural_language_summary(self, tp53_missense: Variant) -> None:
        """The report should have a human-readable summary."""
        report = analyze_variant(tp53_missense)

        assert report.natural_language_summary != ""
        assert "TP53" in report.natural_language_summary or tp53_missense.variant_id in report.natural_language_summary

    def test_variant_without_gene(self) -> None:
        """A variant without a gene symbol should still be processed."""
        variant = Variant(
            chromosome="chr1",
            position=100000,
            reference="A",
            alternate="G",
            variant_type=VariantType.SNV,
        )
        report = analyze_variant(variant)

        assert isinstance(report, TriageReport)
        assert report.classification is not None

    def test_different_variants_get_different_trace_ids(
        self, tp53_missense: Variant, brca1_frameshift: Variant
    ) -> None:
        """Each analysis run should get a unique trace ID."""
        report1 = analyze_variant(tp53_missense)
        report2 = analyze_variant(brca1_frameshift)

        assert report1.trace_id != report2.trace_id


class TestDynamicRouting:
    """Test that conditional edges route correctly."""

    def test_qc_pass_continues_to_annotation(self, tp53_missense: Variant) -> None:
        """When QC passes, the pipeline should continue to annotation."""
        report = analyze_variant(tp53_missense)

        # Check provenance shows annotation happened
        agent_sequence = [p.agent for p in report.provenance]
        assert "annotation_agent" in agent_sequence

    def test_novel_variant_triggers_literature(self, tp53_missense: Variant) -> None:
        """When ClinVar doesn't find the variant, literature search should trigger.

        Since our placeholder annotation returns clinvar.found=False and the
        variant has a gene (TP53), literature should be triggered.
        """
        report = analyze_variant(tp53_missense)

        agent_sequence = [p.agent for p in report.provenance]
        assert "literature_agent" in agent_sequence
