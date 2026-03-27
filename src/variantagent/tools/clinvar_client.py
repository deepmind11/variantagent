"""ClinVar client using NCBI E-utilities.

Queries ClinVar for variant clinical significance using the two-step
esearch → esummary pattern with JSON responses.

Rate limiting: 3 req/sec without API key, 10 req/sec with key.
NCBI does not return rate-limit headers — we self-throttle.
"""

from __future__ import annotations

import asyncio
import logging
from typing import Any

import httpx

from variantagent.config import settings
from variantagent.models.annotation import ClinVarAnnotation
from variantagent.models.variant import Variant

logger = logging.getLogger(__name__)

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Review status → star rating mapping
REVIEW_STATUS_STARS: dict[str, int] = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, single submitter": 1,
    "no assertion criteria provided": 0,
    "no classification provided": 0,
    "no classification for the single variant": 0,
}

# Shared semaphore for NCBI rate limiting (ClinVar + PubMed share the same server)
_ncbi_semaphore: asyncio.Semaphore | None = None


def _get_ncbi_semaphore() -> asyncio.Semaphore:
    """Get or create the NCBI rate-limiting semaphore."""
    global _ncbi_semaphore
    if _ncbi_semaphore is None:
        max_concurrent = 10 if settings.ncbi_api_key else 3
        _ncbi_semaphore = asyncio.Semaphore(max_concurrent)
    return _ncbi_semaphore


def _build_query(variant: Variant) -> str:
    """Build a ClinVar search query from a Variant.

    Tries rsID first (most reliable), then genomic coordinates.
    """
    if variant.rsid:
        return f"{variant.rsid}[rsid]"

    chrom = variant.normalized_chromosome
    return f"{chrom}:{variant.position}:{variant.reference}:{variant.alternate}"


def _base_params() -> dict[str, str]:
    """Common parameters for all NCBI E-utilities requests."""
    params: dict[str, str] = {
        "tool": "variantagent",
        "email": settings.ncbi_email or "variantagent@example.com",
    }
    if settings.ncbi_api_key:
        params["api_key"] = settings.ncbi_api_key
    return params


async def _esearch(client: httpx.AsyncClient, query: str) -> list[str]:
    """Search ClinVar and return matching UIDs."""
    params = {
        **_base_params(),
        "db": "clinvar",
        "term": query,
        "retmode": "json",
        "retmax": "5",
    }

    sem = _get_ncbi_semaphore()
    async with sem:
        response = await client.get(f"{EUTILS_BASE}/esearch.fcgi", params=params)
        response.raise_for_status()

    data = response.json()
    id_list = data.get("esearchresult", {}).get("idlist", [])
    return id_list


async def _esummary(client: httpx.AsyncClient, uids: list[str]) -> dict[str, Any]:
    """Get ClinVar summary data for given UIDs."""
    params = {
        **_base_params(),
        "db": "clinvar",
        "id": ",".join(uids),
        "retmode": "json",
    }

    sem = _get_ncbi_semaphore()
    async with sem:
        response = await client.get(f"{EUTILS_BASE}/esummary.fcgi", params=params)
        response.raise_for_status()

    return response.json()


def _parse_esummary(data: dict[str, Any], uids: list[str]) -> ClinVarAnnotation:
    """Parse esummary JSON response into ClinVarAnnotation."""
    result = data.get("result", {})

    if not uids:
        return ClinVarAnnotation(found=False)

    uid = uids[0]
    record = result.get(uid, {})

    if not record or "error" in record:
        return ClinVarAnnotation(found=False)

    # Extract germline classification
    germline = record.get("germline_classification", {})
    clinical_significance = germline.get("description")
    review_status = germline.get("review_status", "")

    # Extract conditions from trait_set
    conditions: list[str] = []
    trait_sets = germline.get("trait_set", [])
    for trait_set in trait_sets:
        trait_name = trait_set.get("trait_name")
        if trait_name:
            conditions.append(trait_name)

    # Count submitters from supporting_submissions
    submissions = record.get("supporting_submissions", {})
    submitter_count = submissions.get("scv", 0)
    if isinstance(submitter_count, list):
        submitter_count = len(submitter_count)

    # Map review status to stars
    review_stars = REVIEW_STATUS_STARS.get(review_status.lower(), 0) if review_status else 0

    # Last evaluated date
    last_evaluated = germline.get("last_evaluated")

    return ClinVarAnnotation(
        variation_id=str(uid),
        clinical_significance=clinical_significance,
        review_status=review_status,
        review_stars=review_stars,
        conditions=conditions,
        submitter_count=submitter_count if isinstance(submitter_count, int) else 0,
        last_evaluated=last_evaluated,
        found=True,
    )


async def query_clinvar(
    variant: Variant,
    client: httpx.AsyncClient | None = None,
) -> ClinVarAnnotation:
    """Query ClinVar for variant clinical significance.

    Args:
        variant: The variant to look up.
        client: Optional httpx.AsyncClient for connection reuse.

    Returns:
        ClinVarAnnotation with clinical significance, review status, and conditions.
        Returns found=False if the variant is not in ClinVar.
    """
    query = _build_query(variant)
    logger.info("Querying ClinVar: %s", query)

    should_close = client is None
    if client is None:
        client = httpx.AsyncClient(timeout=30.0)

    try:
        uids = await _esearch(client, query)

        if not uids:
            logger.info("Variant not found in ClinVar: %s", variant.variant_id)
            return ClinVarAnnotation(found=False)

        data = await _esummary(client, uids)
        annotation = _parse_esummary(data, uids)

        logger.info(
            "ClinVar result for %s: %s (%s)",
            variant.variant_id,
            annotation.clinical_significance,
            annotation.review_status,
        )
        return annotation

    except httpx.HTTPStatusError as e:
        logger.error("ClinVar HTTP error: %s", e)
        return ClinVarAnnotation(found=False)
    except Exception as e:
        logger.error("ClinVar query failed: %s", e)
        return ClinVarAnnotation(found=False)
    finally:
        if should_close:
            await client.aclose()


async def query_clinvar_safe(
    variant: Variant,
    client: httpx.AsyncClient | None = None,
) -> tuple[ClinVarAnnotation, str | None]:
    """Query ClinVar with error capture — never raises.

    Returns:
        Tuple of (annotation, error_message). error_message is None on success.
    """
    try:
        annotation = await query_clinvar(variant, client)
        return annotation, None
    except Exception as e:
        return ClinVarAnnotation(found=False), f"ClinVar error: {e}"
