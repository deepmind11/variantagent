"""Tests for samtools flagstat parser."""

import pytest

from variantagent.tools.flagstat_parser import parse_flagstat_text


class TestFlagstatParser:
    def test_parse_passing_sample(self, sample_flagstat_text: str) -> None:
        metrics = parse_flagstat_text(sample_flagstat_text)

        assert metrics.total_reads == 10_000_000
        assert metrics.mapped_reads == 9_800_000
        assert metrics.mapping_rate == pytest.approx(0.98, abs=0.01)
        assert metrics.duplicates == 500_000
        assert metrics.duplication_rate == pytest.approx(0.05, abs=0.01)
        assert metrics.properly_paired == 9_500_000
        assert metrics.properly_paired_rate == pytest.approx(0.955, abs=0.01)
        assert metrics.singletons == 100_000
        assert metrics.singleton_rate == pytest.approx(0.01, abs=0.01)

    def test_parse_failing_sample(self, failing_flagstat_text: str) -> None:
        metrics = parse_flagstat_text(failing_flagstat_text)

        assert metrics.total_reads == 5_000_000
        assert metrics.duplicates == 3_000_000
        assert metrics.duplication_rate == pytest.approx(0.60, abs=0.01)
        assert metrics.mapping_rate == pytest.approx(0.80, abs=0.01)
        assert metrics.properly_paired_rate == pytest.approx(0.70, abs=0.01)
        assert metrics.singleton_rate == pytest.approx(0.10, abs=0.01)

    def test_parse_empty_input(self) -> None:
        metrics = parse_flagstat_text("")
        assert metrics.total_reads == 0
        assert metrics.mapping_rate == 0.0
        assert metrics.duplication_rate == 0.0

    def test_parse_zero_reads(self) -> None:
        text = "0 + 0 in total (QC-passed reads + QC-failed reads)\n0 + 0 mapped (0.00% : N/A)"
        metrics = parse_flagstat_text(text)
        assert metrics.total_reads == 0
        assert metrics.mapping_rate == 0.0
