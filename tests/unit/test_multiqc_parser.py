"""Tests for MultiQC parser."""

import json
import tempfile
from pathlib import Path

import pytest

from variantagent.tools.multiqc_parser import (
    _safe_float,
    _safe_int,
    parse_multiqc_data,
    parse_multiqc_json,
)


class TestSafeFloat:
    def test_none_returns_none(self) -> None:
        assert _safe_float(None) is None

    def test_int_converts(self) -> None:
        assert _safe_float(10) == 10.0

    def test_float_passthrough(self) -> None:
        assert _safe_float(3.14) == pytest.approx(3.14)

    def test_string_numeric(self) -> None:
        assert _safe_float("2.5") == pytest.approx(2.5)

    def test_invalid_string_returns_none(self) -> None:
        assert _safe_float("not_a_number") is None

    def test_empty_string_returns_none(self) -> None:
        assert _safe_float("") is None


class TestSafeInt:
    def test_none_returns_none(self) -> None:
        assert _safe_int(None) is None

    def test_int_passthrough(self) -> None:
        assert _safe_int(42) == 42

    def test_float_truncates(self) -> None:
        assert _safe_int(3.9) == 3

    def test_string_numeric(self) -> None:
        assert _safe_int("100") == 100

    def test_invalid_string_returns_none(self) -> None:
        assert _safe_int("abc") is None


class TestParseMultiqcData:
    def test_empty_data_returns_empty_list(self) -> None:
        result = parse_multiqc_data({})
        assert result == []

    def test_list_format_single_block(self) -> None:
        data = {
            "report_general_stats_data": [
                {
                    "sample1": {
                        "total_sequences": 1000000,
                        "percent_gc": 45.0,
                        "mean_coverage": 80.0,
                    }
                }
            ]
        }
        result = parse_multiqc_data(data)
        assert len(result) == 1
        assert result[0].sample_id == "sample1"
        assert result[0].total_sequences == 1000000
        assert result[0].percent_gc == pytest.approx(45.0)
        assert result[0].mean_coverage == pytest.approx(80.0)

    def test_dict_format(self) -> None:
        data = {
            "report_general_stats_data": {
                "sampleA": {
                    "total_sequences": 500000,
                    "percent_duplicates": 10.5,
                },
                "sampleB": {
                    "total_sequences": 750000,
                    "percent_duplicates": 5.0,
                },
            }
        }
        result = parse_multiqc_data(data)
        assert len(result) == 2
        sample_ids = {m.sample_id for m in result}
        assert sample_ids == {"sampleA", "sampleB"}

    def test_general_stats_fallback_key(self) -> None:
        data = {"general_stats": [{"sample1": {"mean_coverage": 50.0}}]}
        result = parse_multiqc_data(data)
        assert len(result) == 1
        assert result[0].mean_coverage == pytest.approx(50.0)

    def test_mosdepth_coverage_field(self) -> None:
        data = {"report_general_stats_data": [{"sample1": {"mosdepth_mean_coverage": 100.0}}]}
        result = parse_multiqc_data(data)
        assert result[0].mean_coverage == pytest.approx(100.0)

    def test_missing_fields_are_none(self) -> None:
        data = {"report_general_stats_data": [{"sample1": {}}]}
        result = parse_multiqc_data(data)
        assert result[0].mean_coverage is None
        assert result[0].total_sequences is None


class TestParseMultiqcJson:
    def test_file_not_found_raises(self) -> None:
        with pytest.raises(FileNotFoundError):
            parse_multiqc_json("/nonexistent/path/multiqc.json")

    def test_valid_json_file(self) -> None:
        data = {
            "report_general_stats_data": [
                {"sample1": {"mean_coverage": 60.0, "total_sequences": 2000000}}
            ]
        }
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(data, f)
            tmp_path = f.name

        result = parse_multiqc_json(tmp_path)
        assert len(result) == 1
        assert result[0].mean_coverage == pytest.approx(60.0)
        Path(tmp_path).unlink()
