"""Tests for the CLI module — testing the pure utility functions."""

import pytest
import typer

from variantagent.cli import _parse_variant_string


class TestParseVariantString:
    def test_format_with_arrow(self) -> None:
        result = _parse_variant_string("chr17:7674220 G>A")
        assert result["chr"] == "chr17"
        assert result["pos"] == "7674220"
        assert result["ref"] == "G"
        assert result["alt"] == "A"

    def test_format_colon_separated(self) -> None:
        result = _parse_variant_string("17:7674220:G:A")
        assert result["chr"] == "17"
        assert result["pos"] == "7674220"
        assert result["ref"] == "G"
        assert result["alt"] == "A"

    def test_format_with_space(self) -> None:
        result = _parse_variant_string("chr17:7674220 G A")
        assert result["chr"] == "chr17"
        assert result["pos"] == "7674220"
        assert result["ref"] == "G"
        assert result["alt"] == "A"

    def test_leading_whitespace_stripped(self) -> None:
        result = _parse_variant_string("  chr1:100000 A>T  ")
        assert result["chr"] == "chr1"
        assert result["ref"] == "A"
        assert result["alt"] == "T"

    def test_multibase_allele(self) -> None:
        result = _parse_variant_string("chr17:43091434 TG>T")
        assert result["ref"] == "TG"
        assert result["alt"] == "T"

    def test_invalid_format_raises(self) -> None:
        with pytest.raises(typer.BadParameter):
            _parse_variant_string("not_a_variant")

    def test_missing_position_raises(self) -> None:
        with pytest.raises(typer.BadParameter):
            _parse_variant_string("chrX G>A")

    def test_case_insensitive(self) -> None:
        result = _parse_variant_string("CHR1:1000:a:t")
        assert result["ref"].upper() == "A"
        assert result["alt"].upper() == "T"
