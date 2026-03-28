"""PubMed client using NCBI E-utilities.

Searches PubMed for variant-specific literature evidence using
progressive broadening: gene + HGVS → gene + protein change → gene + disease.

Shares rate limiting with ClinVar (same NCBI server).
"""

from __future__ import annotations

import logging
from typing import Any

import httpx

from variantagent.config import settings
from variantagent.tools.clinvar_client import _base_params, _get_ncbi_semaphore

logger = logging.getLogger(__name__)

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


class PubMedArticle:
    """A PubMed article summary."""

    def __init__(
        self,
        pmid: str,
        title: str,
        journal: str = "",
        year: str = "",
        authors: list[str] | None = None,
    ) -> None:
        self.pmid = pmid
        self.title = title
        self.journal = journal
        self.year = year
        self.authors = authors or []

    def citation(self) -> str:
        """Format as a short citation."""
        first_author = self.authors[0] if self.authors else "Unknown"
        return f"{first_author} et al. ({self.year}). {self.title}. {self.journal}. PMID: {self.pmid}"


async def _esearch_pubmed(
    client: httpx.AsyncClient,
    query: str,
    max_results: int = 5,
) -> list[str]:
    """Search PubMed and return PMIDs."""
    params = {
        **_base_params(),
        "db": "pubmed",
        "term": query,
        "retmode": "json",
        "retmax": str(max_results),
        "sort": "relevance",
    }

    sem = _get_ncbi_semaphore()
    async with sem:
        response = await client.get(f"{EUTILS_BASE}/esearch.fcgi", params=params)
        response.raise_for_status()

    data = response.json()
    return data.get("esearchresult", {}).get("idlist", [])


async def _esummary_pubmed(
    client: httpx.AsyncClient,
    pmids: list[str],
) -> list[PubMedArticle]:
    """Get article summaries for given PMIDs."""
    if not pmids:
        return []

    params = {
        **_base_params(),
        "db": "pubmed",
        "id": ",".join(pmids),
        "retmode": "json",
    }

    sem = _get_ncbi_semaphore()
    async with sem:
        response = await client.get(f"{EUTILS_BASE}/esummary.fcgi", params=params)
        response.raise_for_status()

    data = response.json()
    result = data.get("result", {})
    articles: list[PubMedArticle] = []

    for pmid in pmids:
        record = result.get(pmid, {})
        if not record or "error" in record:
            continue

        # Parse authors
        authors = []
        for author in record.get("authors", []):
            name = author.get("name", "")
            if name:
                authors.append(name)

        articles.append(
            PubMedArticle(
                pmid=pmid,
                title=record.get("title", ""),
                journal=record.get("fulljournalname", record.get("source", "")),
                year=record.get("pubdate", "")[:4],
                authors=authors,
            )
        )

    return articles


def _build_search_queries(
    gene: str | None,
    hgvs_p: str | None = None,
    variant_id: str | None = None,
) -> list[str]:
    """Build progressively broader PubMed search queries.

    Strategy: specific → broad
    1. Gene + protein change (most specific)
    2. Gene + "variant" + "pathogenic" (broader)
    3. Gene + "ACMG" (ACMG-specific literature)
    """
    queries: list[str] = []

    if gene and hgvs_p:
        # Most specific: gene + protein change
        protein_change = hgvs_p.replace("p.", "")
        queries.append(f"{gene}[Gene] AND {protein_change}")

    if gene:
        queries.append(f"{gene}[Gene] AND (variant classification OR pathogenic OR ACMG)")
        queries.append(f"{gene}[Gene] AND variant interpretation")

    return queries


async def search_pubmed(
    gene: str | None,
    hgvs_p: str | None = None,
    variant_id: str | None = None,
    max_results: int = 5,
    client: httpx.AsyncClient | None = None,
) -> list[PubMedArticle]:
    """Search PubMed for variant-relevant literature.

    Uses progressive broadening: tries the most specific query first,
    falls back to broader queries if insufficient results.

    Args:
        gene: Gene symbol (e.g., "TP53").
        hgvs_p: HGVS protein notation (e.g., "p.R175H").
        variant_id: Variant ID for logging.
        max_results: Maximum number of articles to return.
        client: Optional httpx.AsyncClient for connection reuse.

    Returns:
        List of PubMedArticle objects, sorted by relevance.
    """
    if not gene:
        logger.info("No gene specified — skipping PubMed search")
        return []

    queries = _build_search_queries(gene, hgvs_p, variant_id)
    logger.info("PubMed search queries for %s: %s", gene, queries)

    should_close = client is None
    if client is None:
        client = httpx.AsyncClient(timeout=30.0)

    try:
        all_pmids: list[str] = []
        seen: set[str] = set()

        for query in queries:
            if len(all_pmids) >= max_results:
                break

            pmids = await _esearch_pubmed(client, query, max_results=max_results)
            for pmid in pmids:
                if pmid not in seen and len(all_pmids) < max_results:
                    all_pmids.append(pmid)
                    seen.add(pmid)

        if not all_pmids:
            logger.info("No PubMed results for %s", gene)
            return []

        articles = await _esummary_pubmed(client, all_pmids)
        logger.info("Found %d PubMed articles for %s", len(articles), gene)
        return articles

    except Exception as e:
        logger.error("PubMed search failed: %s", e)
        return []
    finally:
        if should_close:
            await client.aclose()


async def search_pubmed_safe(
    gene: str | None,
    hgvs_p: str | None = None,
    variant_id: str | None = None,
    max_results: int = 5,
    client: httpx.AsyncClient | None = None,
) -> tuple[list[PubMedArticle], str | None]:
    """Search PubMed with error capture — never raises."""
    try:
        articles = await search_pubmed(gene, hgvs_p, variant_id, max_results, client)
        return articles, None
    except Exception as e:
        return [], f"PubMed error: {e}"
