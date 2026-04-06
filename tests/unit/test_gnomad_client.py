"""Tests for gnomAD client — parsing logic only (no live API calls)."""

import pytest

from variantagent.models.variant import Variant
from variantagent.tools.gnomad_client import (
    GNOMAD_API,
    VARIANT_QUERY,
    _build_variant_id,
    _compute_af,
    _parse_gnomad_response,
)


class TestBuildVariantId:
    def test_standard_variant(self) -> None:
        v = Variant(chromosome="chr17", position=7674220, reference="G", alternate="A")
        assert _build_variant_id(v) == "17-7674220-G-A"

    def test_strips_chr_prefix(self) -> None:
        v = Variant(chromosome="chr1", position=100, reference="A", alternate="G")
        vid = _build_variant_id(v)
        assert vid == "1-100-A-G"
        assert "chr" not in vid


class TestComputeAf:
    def test_normal_af(self) -> None:
        assert _compute_af(10, 1000) == pytest.approx(0.01)

    def test_zero_an(self) -> None:
        assert _compute_af(10, 0) is None

    def test_none_values(self) -> None:
        assert _compute_af(None, 1000) is None
        assert _compute_af(10, None) is None


class TestParseGnomadResponse:
    def test_parse_exome_data(self) -> None:
        data = {
            "data": {
                "variant": {
                    "exome": {
                        "ac": 50,
                        "an": 100000,
                        "ac_hom": 2,
                        "populations": [
                            {"id": "afr", "ac": 10, "an": 20000},
                            {"id": "nfe", "ac": 30, "an": 60000},
                            {"id": "eas", "ac": 5, "an": 10000},
                        ],
                    },
                    "genome": None,
                }
            }
        }
        result = _parse_gnomad_response(data)
        assert result.found is True
        assert result.overall_af == pytest.approx(0.0005)
        assert result.afr_af == pytest.approx(0.0005)
        assert result.nfe_af == pytest.approx(0.0005)
        assert result.eas_af == pytest.approx(0.0005)
        assert result.homozygote_count == 2
        assert result.allele_count == 50
        assert result.allele_number == 100000

    def test_parse_combined_exome_genome(self) -> None:
        data = {
            "data": {
                "variant": {
                    "exome": {"ac": 30, "an": 50000, "ac_hom": 1, "populations": []},
                    "genome": {"ac": 20, "an": 30000, "ac_hom": 0, "populations": []},
                }
            }
        }
        result = _parse_gnomad_response(data)
        assert result.found is True
        assert result.allele_count == 50
        assert result.allele_number == 80000
        assert result.homozygote_count == 1
        assert result.overall_af == pytest.approx(50 / 80000)

    def test_parse_not_found(self) -> None:
        data = {"data": {"variant": None}}
        result = _parse_gnomad_response(data)
        assert result.found is False

    def test_parse_both_null(self) -> None:
        data = {"data": {"variant": {"exome": None, "genome": None}}}
        result = _parse_gnomad_response(data)
        assert result.found is False

    def test_parse_genome_only(self) -> None:
        """Only genome data (no exome)."""
        data = {
            "data": {
                "variant": {
                    "exome": None,
                    "genome": {
                        "ac": 5,
                        "an": 10000,
                        "ac_hom": 0,
                        "populations": [{"id": "nfe", "ac": 4, "an": 8000}],
                    },
                }
            }
        }
        result = _parse_gnomad_response(data)
        assert result.found is True
        assert result.allele_count == 5
        assert result.allele_number == 10000
        assert result.nfe_af == pytest.approx(4 / 8000)

    def test_parse_population_merged_from_both(self) -> None:
        """Population data merged from exome + genome."""
        data = {
            "data": {
                "variant": {
                    "exome": {
                        "ac": 10,
                        "an": 50000,
                        "ac_hom": 0,
                        "populations": [{"id": "afr", "ac": 5, "an": 10000}],
                    },
                    "genome": {
                        "ac": 5,
                        "an": 30000,
                        "ac_hom": 0,
                        "populations": [{"id": "afr", "ac": 3, "an": 6000}],
                    },
                }
            }
        }
        result = _parse_gnomad_response(data)
        assert result.found is True
        # afr merged: ac=8, an=16000 → af=0.0005
        assert result.afr_af == pytest.approx(8 / 16000)


class TestConstants:
    def test_gnomad_api_url(self) -> None:
        assert "gnomad" in GNOMAD_API
        assert GNOMAD_API.startswith("https://")

    def test_variant_query_is_graphql(self) -> None:
        assert "query" in VARIANT_QUERY.lower()
        assert "variantId" in VARIANT_QUERY
