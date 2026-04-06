"""Tests for ClinVar client — parsing logic only (no live API calls)."""

import asyncio

from variantagent.models.variant import Variant
from variantagent.tools.clinvar_client import (
    EUTILS_BASE,
    REVIEW_STATUS_STARS,
    _base_params,
    _build_query,
    _get_ncbi_semaphore,
    _parse_esummary,
)


class TestBuildQuery:
    def test_query_by_rsid(self) -> None:
        v = Variant(
            chromosome="chr17",
            position=7674220,
            reference="G",
            alternate="A",
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

    def test_parse_record_without_conditions(self) -> None:
        data = {
            "result": {
                "uids": ["11111"],
                "11111": {
                    "uid": "11111",
                    "germline_classification": {
                        "description": "Likely pathogenic",
                        "review_status": "criteria provided, single submitter",
                        "trait_set": [],
                    },
                    "supporting_submissions": {"scv": ["SCV001"]},
                },
            }
        }
        result = _parse_esummary(data, ["11111"])
        assert result.found is True
        assert result.conditions == []

    def test_parse_expert_panel_review_stars(self) -> None:
        data = {
            "result": {
                "uids": ["22222"],
                "22222": {
                    "uid": "22222",
                    "germline_classification": {
                        "description": "Pathogenic",
                        "review_status": "reviewed by expert panel",
                        "trait_set": [],
                    },
                    "supporting_submissions": {"scv": []},
                },
            }
        }
        result = _parse_esummary(data, ["22222"])
        assert result.review_stars == 3


class TestBaseParams:
    def test_returns_tool_and_email(self) -> None:
        params = _base_params()
        assert params["tool"] == "variantagent"
        assert "email" in params

    def test_returns_dict_of_strings(self) -> None:
        params = _base_params()
        assert isinstance(params, dict)
        for k, v in params.items():
            assert isinstance(k, str)
            assert isinstance(v, str)


class TestGetNcbiSemaphore:
    def test_returns_semaphore(self) -> None:
        import variantagent.tools.clinvar_client as cc

        cc._ncbi_semaphore = None  # reset
        sem = _get_ncbi_semaphore()
        assert isinstance(sem, asyncio.Semaphore)

    def test_returns_same_instance(self) -> None:
        sem1 = _get_ncbi_semaphore()
        sem2 = _get_ncbi_semaphore()
        assert sem1 is sem2


class TestReviewStatusStars:
    def test_practice_guideline_is_4(self) -> None:
        assert REVIEW_STATUS_STARS["practice guideline"] == 4

    def test_expert_panel_is_3(self) -> None:
        assert REVIEW_STATUS_STARS["reviewed by expert panel"] == 3

    def test_no_criteria_is_0(self) -> None:
        assert REVIEW_STATUS_STARS["no assertion criteria provided"] == 0


class TestEutilsBase:
    def test_eutils_base_is_ncbi_url(self) -> None:
        assert "ncbi.nlm.nih.gov" in EUTILS_BASE
        assert EUTILS_BASE.startswith("https://")
