"""Async tests for client wrapper functions using mocks."""

from unittest.mock import AsyncMock, MagicMock, patch

import pytest

from variantagent.models.annotation import ClinVarAnnotation, EnsemblVEPAnnotation, GnomADFrequency
from variantagent.models.variant import Variant


@pytest.fixture
def variant() -> Variant:
    return Variant(
        chromosome="chr17",
        position=7674220,
        reference="G",
        alternate="A",
        gene="TP53",
        rsid="rs28934578",
    )


class TestClinVarSafe:
    async def test_safe_returns_error_message_on_exception(self, variant: Variant) -> None:
        from variantagent.tools.clinvar_client import query_clinvar_safe

        with patch(
            "variantagent.tools.clinvar_client.query_clinvar",
            new_callable=AsyncMock,
            side_effect=RuntimeError("network error"),
        ):
            annotation, error = await query_clinvar_safe(variant)
            assert annotation.found is False
            assert error is not None
            assert "ClinVar error" in error

    async def test_safe_returns_none_error_on_success(self, variant: Variant) -> None:
        from variantagent.tools.clinvar_client import query_clinvar_safe

        mock_result = ClinVarAnnotation(found=True, clinical_significance="Pathogenic")
        with patch(
            "variantagent.tools.clinvar_client.query_clinvar",
            new_callable=AsyncMock,
            return_value=mock_result,
        ):
            annotation, error = await query_clinvar_safe(variant)
            assert annotation.found is True
            assert error is None


class TestGnomADSafe:
    async def test_safe_returns_error_message_on_exception(self, variant: Variant) -> None:
        from variantagent.tools.gnomad_client import query_gnomad_safe

        with patch(
            "variantagent.tools.gnomad_client.query_gnomad",
            new_callable=AsyncMock,
            side_effect=RuntimeError("timeout"),
        ):
            freq, error = await query_gnomad_safe(variant)
            assert freq.found is False
            assert error is not None
            assert "gnomAD error" in error

    async def test_safe_returns_none_error_on_success(self, variant: Variant) -> None:
        from variantagent.tools.gnomad_client import query_gnomad_safe

        mock_result = GnomADFrequency(found=True, overall_af=0.001)
        with patch(
            "variantagent.tools.gnomad_client.query_gnomad",
            new_callable=AsyncMock,
            return_value=mock_result,
        ):
            freq, error = await query_gnomad_safe(variant)
            assert freq.found is True
            assert error is None


class TestVepSafe:
    async def test_safe_returns_error_message_on_exception(self, variant: Variant) -> None:
        from variantagent.tools.ensembl_client import query_vep_safe

        with patch(
            "variantagent.tools.ensembl_client.query_vep",
            new_callable=AsyncMock,
            side_effect=RuntimeError("connection refused"),
        ):
            annotation, error = await query_vep_safe(variant)
            assert annotation.found is False
            assert error is not None
            assert "VEP error" in error

    async def test_safe_returns_none_error_on_success(self, variant: Variant) -> None:
        from variantagent.tools.ensembl_client import query_vep_safe

        mock_result = EnsemblVEPAnnotation(found=True, consequence_type="missense_variant")
        with patch(
            "variantagent.tools.ensembl_client.query_vep",
            new_callable=AsyncMock,
            return_value=mock_result,
        ):
            annotation, error = await query_vep_safe(variant)
            assert annotation.found is True
            assert error is None


class TestClinVarQueryMocked:
    async def test_query_clinvar_not_found_when_no_uids(self, variant: Variant) -> None:
        from variantagent.tools.clinvar_client import query_clinvar

        with patch(
            "variantagent.tools.clinvar_client._esearch",
            new_callable=AsyncMock,
            return_value=[],
        ):
            result = await query_clinvar(variant)
            assert result.found is False

    async def test_query_clinvar_http_error_returns_not_found(self, variant: Variant) -> None:
        import httpx

        from variantagent.tools.clinvar_client import query_clinvar

        mock_response = MagicMock()
        mock_response.status_code = 429
        http_error = httpx.HTTPStatusError(
            "rate limited", request=MagicMock(), response=mock_response
        )
        with patch(
            "variantagent.tools.clinvar_client._esearch",
            new_callable=AsyncMock,
            side_effect=http_error,
        ):
            result = await query_clinvar(variant)
            assert result.found is False

    async def test_query_clinvar_found(self, variant: Variant) -> None:
        from variantagent.tools.clinvar_client import query_clinvar

        mock_annotation = ClinVarAnnotation(
            found=True, clinical_significance="Pathogenic", review_stars=2
        )
        with (
            patch(
                "variantagent.tools.clinvar_client._esearch",
                new_callable=AsyncMock,
                return_value=["12345"],
            ),
            patch(
                "variantagent.tools.clinvar_client._esummary",
                new_callable=AsyncMock,
                return_value={"result": {}},
            ),
            patch(
                "variantagent.tools.clinvar_client._parse_esummary",
                return_value=mock_annotation,
            ),
        ):
            result = await query_clinvar(variant)
            assert result.found is True
            assert result.clinical_significance == "Pathogenic"


class TestGnomADQueryMocked:
    async def test_query_gnomad_found(self, variant: Variant) -> None:
        from variantagent.tools.gnomad_client import query_gnomad

        response_data = {
            "data": {
                "variant": {
                    "exome": {"ac": 5, "an": 10000, "ac_hom": 0, "populations": []},
                    "genome": None,
                }
            }
        }
        mock_response = MagicMock()
        mock_response.raise_for_status = MagicMock()
        mock_response.json = MagicMock(return_value=response_data)

        import httpx

        mock_client = AsyncMock(spec=httpx.AsyncClient)
        mock_client.post = AsyncMock(return_value=mock_response)

        with patch("variantagent.tools.gnomad_client.asyncio.sleep", new_callable=AsyncMock):
            result = await query_gnomad(variant, client=mock_client)
            assert result.found is True
            assert result.allele_count == 5

    async def test_query_gnomad_graphql_errors(self, variant: Variant) -> None:
        from variantagent.tools.gnomad_client import query_gnomad

        mock_response = MagicMock()
        mock_response.raise_for_status = MagicMock()
        mock_response.json = MagicMock(return_value={"errors": [{"message": "variant not found"}]})

        import httpx

        mock_client = AsyncMock(spec=httpx.AsyncClient)
        mock_client.post = AsyncMock(return_value=mock_response)

        with patch("variantagent.tools.gnomad_client.asyncio.sleep", new_callable=AsyncMock):
            result = await query_gnomad(variant, client=mock_client)
            assert result.found is False

    async def test_query_gnomad_not_found(self, variant: Variant) -> None:
        from variantagent.tools.gnomad_client import query_gnomad

        mock_response = MagicMock()
        mock_response.raise_for_status = MagicMock()
        mock_response.json = MagicMock(return_value={"data": {"variant": None}})

        import httpx

        mock_client = AsyncMock(spec=httpx.AsyncClient)
        mock_client.post = AsyncMock(return_value=mock_response)

        with patch("variantagent.tools.gnomad_client.asyncio.sleep", new_callable=AsyncMock):
            result = await query_gnomad(variant, client=mock_client)
            assert result.found is False


class TestVepQueryMocked:
    async def test_query_vep_found(self, variant: Variant) -> None:
        from variantagent.tools.ensembl_client import query_vep

        vep_data = [
            {
                "most_severe_consequence": "missense_variant",
                "transcript_consequences": [
                    {
                        "consequence_terms": ["missense_variant"],
                        "impact": "MODERATE",
                        "gene_symbol": "TP53",
                    }
                ],
            }
        ]
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.raise_for_status = MagicMock()
        mock_response.json = MagicMock(return_value=vep_data)

        import httpx

        mock_client = AsyncMock(spec=httpx.AsyncClient)
        mock_client.get = AsyncMock(return_value=mock_response)

        result = await query_vep(variant, client=mock_client)
        assert result.found is True
        assert result.consequence_type == "missense_variant"

    async def test_query_vep_rate_limited(self, variant: Variant) -> None:
        from variantagent.tools.ensembl_client import query_vep

        mock_response = MagicMock()
        mock_response.status_code = 429
        mock_response.headers = {"Retry-After": "5"}

        import httpx

        mock_client = AsyncMock(spec=httpx.AsyncClient)
        mock_client.get = AsyncMock(return_value=mock_response)

        result = await query_vep(variant, client=mock_client)
        assert result.found is False

    async def test_query_vep_http_error(self, variant: Variant) -> None:
        import httpx

        from variantagent.tools.ensembl_client import query_vep

        mock_response = MagicMock()
        mock_response.status_code = 500
        http_error = httpx.HTTPStatusError(
            "server error", request=MagicMock(), response=mock_response
        )

        mock_client = AsyncMock(spec=httpx.AsyncClient)
        mock_client.get = AsyncMock(side_effect=http_error)

        result = await query_vep(variant, client=mock_client)
        assert result.found is False


class TestPubMedSearch:
    async def test_search_no_gene_returns_empty(self) -> None:
        from variantagent.tools.pubmed_client import search_pubmed

        result = await search_pubmed(None)
        assert result == []

    async def test_search_found_articles(self) -> None:
        from variantagent.tools.pubmed_client import search_pubmed

        mock_articles = [
            MagicMock(pmid="111", title="TP53 study"),
        ]
        with (
            patch(
                "variantagent.tools.pubmed_client._esearch_pubmed",
                new_callable=AsyncMock,
                return_value=["111"],
            ),
            patch(
                "variantagent.tools.pubmed_client._esummary_pubmed",
                new_callable=AsyncMock,
                return_value=mock_articles,
            ),
        ):
            import httpx

            async with httpx.AsyncClient() as client:
                result = await search_pubmed("TP53", client=client)
                assert len(result) == 1

    async def test_search_no_pmids_returns_empty(self) -> None:
        from variantagent.tools.pubmed_client import search_pubmed

        with patch(
            "variantagent.tools.pubmed_client._esearch_pubmed",
            new_callable=AsyncMock,
            return_value=[],
        ):
            import httpx

            async with httpx.AsyncClient() as client:
                result = await search_pubmed("BRCA1", client=client)
                assert result == []

    async def test_search_exception_returns_empty(self) -> None:
        from variantagent.tools.pubmed_client import search_pubmed

        with patch(
            "variantagent.tools.pubmed_client._esearch_pubmed",
            new_callable=AsyncMock,
            side_effect=RuntimeError("network error"),
        ):
            import httpx

            async with httpx.AsyncClient() as client:
                result = await search_pubmed("TP53", client=client)
                assert result == []

    async def test_safe_search_returns_error_string(self) -> None:
        from variantagent.tools.pubmed_client import search_pubmed_safe

        with patch(
            "variantagent.tools.pubmed_client.search_pubmed",
            new_callable=AsyncMock,
            side_effect=RuntimeError("pubmed down"),
        ):
            articles, error = await search_pubmed_safe("TP53")
            assert articles == []
            assert error is not None

    async def test_esummary_pubmed_empty_pmids(self) -> None:
        """_esummary_pubmed returns empty list for empty PMIDs."""
        import httpx

        from variantagent.tools.pubmed_client import _esummary_pubmed

        async with httpx.AsyncClient() as client:
            result = await _esummary_pubmed(client, [])
            assert result == []

    async def test_esummary_pubmed_parses_articles(self) -> None:
        import httpx

        from variantagent.tools.pubmed_client import _esummary_pubmed

        response_data = {
            "result": {
                "111": {
                    "title": "BRCA1 variant analysis",
                    "fulljournalname": "Nature Genetics",
                    "pubdate": "2023 Jan",
                    "authors": [{"name": "Smith AB"}, {"name": "Jones CD"}],
                }
            }
        }
        mock_response = MagicMock()
        mock_response.raise_for_status = MagicMock()
        mock_response.json = MagicMock(return_value=response_data)

        mock_client = AsyncMock(spec=httpx.AsyncClient)
        mock_client.get = AsyncMock(return_value=mock_response)

        with patch("variantagent.tools.pubmed_client._get_ncbi_semaphore") as mock_sem:
            mock_sem.return_value = MagicMock()
            mock_sem.return_value.__aenter__ = AsyncMock(return_value=None)
            mock_sem.return_value.__aexit__ = AsyncMock(return_value=None)
            result = await _esummary_pubmed(mock_client, ["111"])

        assert len(result) == 1
        assert result[0].title == "BRCA1 variant analysis"
        assert result[0].year == "2023"
        assert "Smith AB" in result[0].authors

    async def test_esearch_pubmed_returns_pmids(self) -> None:
        import httpx

        from variantagent.tools.pubmed_client import _esearch_pubmed

        response_data = {"esearchresult": {"idlist": ["111", "222", "333"]}}
        mock_response = MagicMock()
        mock_response.raise_for_status = MagicMock()
        mock_response.json = MagicMock(return_value=response_data)

        mock_client = AsyncMock(spec=httpx.AsyncClient)
        mock_client.get = AsyncMock(return_value=mock_response)

        with patch("variantagent.tools.pubmed_client._get_ncbi_semaphore") as mock_sem:
            mock_sem.return_value = MagicMock()
            mock_sem.return_value.__aenter__ = AsyncMock(return_value=None)
            mock_sem.return_value.__aexit__ = AsyncMock(return_value=None)
            result = await _esearch_pubmed(mock_client, "TP53[Gene] AND R175H")

        assert result == ["111", "222", "333"]


class TestEsearchEsummaryMocked:
    async def test_esearch_returns_ids(self, variant: Variant) -> None:
        import httpx

        from variantagent.tools.clinvar_client import _esearch

        mock_response = MagicMock()
        mock_response.raise_for_status = MagicMock()
        mock_response.json = MagicMock(return_value={"esearchresult": {"idlist": ["111", "222"]}})

        mock_client = AsyncMock(spec=httpx.AsyncClient)
        mock_client.get = AsyncMock(return_value=mock_response)

        with patch("variantagent.tools.clinvar_client._get_ncbi_semaphore") as mock_sem:
            mock_sem.return_value = MagicMock()
            mock_sem.return_value.__aenter__ = AsyncMock(return_value=None)
            mock_sem.return_value.__aexit__ = AsyncMock(return_value=None)
            result = await _esearch(mock_client, "rs12345[rsid]")
            assert result == ["111", "222"]

    async def test_esummary_returns_dict(self, variant: Variant) -> None:
        import httpx

        from variantagent.tools.clinvar_client import _esummary

        expected_data = {"result": {"uids": ["111"], "111": {"uid": "111"}}}
        mock_response = MagicMock()
        mock_response.raise_for_status = MagicMock()
        mock_response.json = MagicMock(return_value=expected_data)

        mock_client = AsyncMock(spec=httpx.AsyncClient)
        mock_client.get = AsyncMock(return_value=mock_response)

        with patch("variantagent.tools.clinvar_client._get_ncbi_semaphore") as mock_sem:
            mock_sem.return_value = MagicMock()
            mock_sem.return_value.__aenter__ = AsyncMock(return_value=None)
            mock_sem.return_value.__aexit__ = AsyncMock(return_value=None)
            result = await _esummary(mock_client, ["111"])
            assert result == expected_data
