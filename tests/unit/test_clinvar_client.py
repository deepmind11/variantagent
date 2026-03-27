"""Tests for ClinVar client — parsing logic only (no live API calls)."""

from variantagent.models.annotation import ClinVarAnnotation
from variantagent.tools.clinvar_client import _build_query, _parse_esummary
from variantagent.models.variant import Variant, VariantType


class TestBuildQuery:
    def test_query_by_rsid(self) -> None:
        v = Variant(
            chromosome="chr17", position=7674220, reference="G", alternate="A",
            rsid="rs28934578",
        )
        assert _build_query(v) == "rs28934578[rsid]"

    def test_query_by_coordinates(self) -> None:
        v = Variant(chromosome="chr17", position=7674220, reference="G", alternate="A")
        assert _build_query(v) == "17:7674220:G:A"

    def test_query_strips_chr_prefix(self) -> None:
        v = Variant(chromosome="chr1", position=100, reference="A", alternate="G")
        query = _build_query(v)
        assert not query.startswith("chr")


class TestParseEsummary:
    def test_parse_pathogenic_variant(self) -> None:
        data = {
            "result": {
                "uids": ["65533"],
                "65533": {
                    "uid": "65533",
                    "accession": "VCV000065533.11",
                    "germline_classification": {
                        "description": "Pathogenic",
                        "last_evaluated": "2025/10/18",
                        "review_status": "criteria provided, multiple submitters, no conflicts",
                        "trait_set": [
                            {"trait_name": "Hereditary breast cancer"},
                            {"trait_name": "Li-Fraumeni syndrome"},
                        ],
                    },
                    "supporting_submissions": {"scv": ["SCV001", "SCV002", "SCV003"]},
                },
            }
        }
        result = _parse_esummary(data, ["65533"])
        assert result.found is True
        assert result.clinical_significance == "Pathogenic"
        assert result.review_stars == 2
        assert "Hereditary breast cancer" in result.conditions
        assert result.submitter_count == 3

    def test_parse_not_found(self) -> None:
        result = _parse_esummary({"result": {}}, [])
        assert result.found is False

    def test_parse_error_record(self) -> None:
        data = {"result": {"uids": ["999"], "999": {"error": "not found"}}}
        result = _parse_esummary(data, ["999"])
        assert result.found is False

    def test_parse_vus(self) -> None:
        data = {
            "result": {
                "uids": ["12345"],
                "12345": {
                    "uid": "12345",
                    "germline_classification": {
                        "description": "Uncertain significance",
                        "review_status": "criteria provided, single submitter",
                        "trait_set": [{"trait_name": "Cardiomyopathy"}],
                    },
                    "supporting_submissions": {"scv": ["SCV100"]},
                },
            }
        }
        result = _parse_esummary(data, ["12345"])
        assert result.found is True
        assert result.clinical_significance == "Uncertain significance"
        assert result.review_stars == 1
