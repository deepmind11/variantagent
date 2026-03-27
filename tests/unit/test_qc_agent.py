"""Tests for the QC Agent assessment logic."""

from variantagent.agents.qc_agent import assess_flagstat, assess_multiqc, run_qc_assessment
from variantagent.models.qc_metrics import FlagstatMetrics, MultiQCMetrics, QCStatus


class TestAssessFlagstat:
    def test_passing_sample_no_issues(self) -> None:
        flagstat = FlagstatMetrics(
            total_reads=10_000_000,
            mapped_reads=9_800_000,
            mapping_rate=0.98,
            duplicates=500_000,
            duplication_rate=0.05,
            paired_reads=9_950_000,
            properly_paired=9_500_000,
            properly_paired_rate=0.955,
            singletons=100_000,
            singleton_rate=0.01,
        )
        issues = assess_flagstat(flagstat)
        assert len(issues) == 0

    def test_high_duplication_fail(self) -> None:
        flagstat = FlagstatMetrics(
            total_reads=5_000_000,
            mapped_reads=4_500_000,
            mapping_rate=0.96,
            duplicates=3_000_000,
            duplication_rate=0.60,
            paired_reads=5_000_000,
            properly_paired=4_500_000,
            properly_paired_rate=0.90,
            singletons=50_000,
            singleton_rate=0.01,
        )
        issues = assess_flagstat(flagstat)
        dup_issues = [i for i in issues if i.metric == "duplication_rate"]
        assert len(dup_issues) == 1
        assert dup_issues[0].severity == QCStatus.FAIL
        assert dup_issues[0].observed_value == 0.60

    def test_low_mapping_rate_fail(self) -> None:
        flagstat = FlagstatMetrics(
            total_reads=5_000_000,
            mapped_reads=4_000_000,
            mapping_rate=0.80,
            duplicates=100_000,
            duplication_rate=0.02,
            paired_reads=5_000_000,
            properly_paired=3_500_000,
            properly_paired_rate=0.70,
            singletons=500_000,
            singleton_rate=0.10,
        )
        issues = assess_flagstat(flagstat)
        map_issues = [i for i in issues if i.metric == "mapping_rate"]
        assert len(map_issues) == 1
        assert map_issues[0].severity == QCStatus.FAIL

    def test_warn_level_duplication(self) -> None:
        flagstat = FlagstatMetrics(
            total_reads=10_000_000,
            mapped_reads=9_800_000,
            mapping_rate=0.98,
            duplicates=3_500_000,
            duplication_rate=0.35,
            paired_reads=10_000_000,
            properly_paired=9_500_000,
            properly_paired_rate=0.95,
            singletons=100_000,
            singleton_rate=0.01,
        )
        issues = assess_flagstat(flagstat)
        dup_issues = [i for i in issues if i.metric == "duplication_rate"]
        assert len(dup_issues) == 1
        assert dup_issues[0].severity == QCStatus.WARN

    def test_multiple_failures(self) -> None:
        """A really bad sample should flag multiple issues."""
        flagstat = FlagstatMetrics(
            total_reads=5_000_000,
            mapped_reads=4_000_000,
            mapping_rate=0.80,
            duplicates=3_000_000,
            duplication_rate=0.60,
            paired_reads=5_000_000,
            properly_paired=3_500_000,
            properly_paired_rate=0.70,
            singletons=600_000,
            singleton_rate=0.12,
        )
        issues = assess_flagstat(flagstat)
        assert len(issues) >= 3  # mapping, duplication, properly_paired, singleton
        fail_count = sum(1 for i in issues if i.severity == QCStatus.FAIL)
        assert fail_count >= 3


class TestAssessMultiqc:
    def test_low_coverage_fail(self) -> None:
        multiqc = MultiQCMetrics(sample_id="S001", mean_coverage=15.0)
        issues = assess_multiqc(multiqc)
        assert len(issues) == 1
        assert issues[0].metric == "mean_coverage"
        assert issues[0].severity == QCStatus.FAIL

    def test_adequate_coverage_no_issues(self) -> None:
        multiqc = MultiQCMetrics(sample_id="S001", mean_coverage=200.0)
        issues = assess_multiqc(multiqc)
        assert len(issues) == 0

    def test_adapter_contamination_warn(self) -> None:
        multiqc = MultiQCMetrics(sample_id="S001", mean_coverage=200.0, percent_adapter=8.0)
        issues = assess_multiqc(multiqc)
        adapter_issues = [i for i in issues if i.metric == "percent_adapter"]
        assert len(adapter_issues) == 1
        assert adapter_issues[0].severity == QCStatus.WARN


class TestRunQCAssessment:
    def test_overall_pass(self) -> None:
        flagstat = FlagstatMetrics(
            total_reads=10_000_000, mapped_reads=9_800_000, mapping_rate=0.98,
            duplicates=500_000, duplication_rate=0.05, paired_reads=10_000_000,
            properly_paired=9_500_000, properly_paired_rate=0.95,
            singletons=100_000, singleton_rate=0.01,
        )
        assessment = run_qc_assessment("S001", flagstat=flagstat, variant_region_coverage=150.0)
        assert assessment.overall_status == QCStatus.PASS
        assert assessment.reliable_for_interpretation is True
        assert len(assessment.issues) == 0

    def test_overall_fail_from_flagstat(self) -> None:
        flagstat = FlagstatMetrics(
            total_reads=5_000_000, mapped_reads=4_000_000, mapping_rate=0.80,
            duplicates=3_000_000, duplication_rate=0.60, paired_reads=5_000_000,
            properly_paired=3_500_000, properly_paired_rate=0.70,
            singletons=600_000, singleton_rate=0.12,
        )
        assessment = run_qc_assessment("S001", flagstat=flagstat)
        assert assessment.overall_status == QCStatus.FAIL
        assert assessment.reliable_for_interpretation is False

    def test_low_variant_coverage_unreliable(self) -> None:
        flagstat = FlagstatMetrics(
            total_reads=10_000_000, mapped_reads=9_800_000, mapping_rate=0.98,
            duplicates=500_000, duplication_rate=0.05, paired_reads=10_000_000,
            properly_paired=9_500_000, properly_paired_rate=0.95,
            singletons=100_000, singleton_rate=0.01,
        )
        assessment = run_qc_assessment("S001", flagstat=flagstat, variant_region_coverage=5.0)
        assert assessment.overall_status == QCStatus.FAIL
        assert assessment.reliable_for_interpretation is False

    def test_no_data_passes(self) -> None:
        """No QC data available — should pass (can't fail what you can't measure)."""
        assessment = run_qc_assessment("S001")
        assert assessment.overall_status == QCStatus.PASS
        assert assessment.reliable_for_interpretation is True

    def test_reasoning_populated(self) -> None:
        assessment = run_qc_assessment("S001")
        assert assessment.reasoning != ""
        assert "reliable" in assessment.reasoning.lower()
