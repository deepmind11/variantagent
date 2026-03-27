"""Ensembl VEP REST API client.

Queries the Variant Effect Predictor for consequence prediction,
SIFT/PolyPhen scores, protein domain information, and impact assessment.

Rate limit: 55,000 requests/hour (~15/sec). Returns Retry-After headers on 429.
"""

from __future__ import annotations

import logging

import httpx

from variantagent.models.annotation import EnsemblVEPAnnotation
from variantagent.models.variant import Variant

logger = logging.getLogger(__name__)

VEP_BASE = "https://rest.ensembl.org"


def _build_vep_url(variant: Variant) -> str:
    """Build the VEP REST API URL for a variant.

    Uses the region-based endpoint: /vep/human/region/{chr}:{pos}:{pos}/{alt}
    """
    chrom = variant.normalized_chromosome
    return (
        f"{VEP_BASE}/vep/human/region/"
        f"{chrom}:{variant.position}:{variant.position}/{variant.alternate}"
    )


def _parse_vep_response(data: list[dict]) -> EnsemblVEPAnnotation:
    """Parse VEP JSON response into EnsemblVEPAnnotation.

    Uses the first result with `pick=1` to get the most severe consequence.
    """
    if not data:
        return EnsemblVEPAnnotation(found=False)

    record = data[0]
    consequences = record.get("transcript_consequences", [])

    if not consequences:
        # Try intergenic consequences
        intergenic = record.get("intergenic_consequences", [])
        if intergenic:
            return EnsemblVEPAnnotation(
                consequence_type=record.get("most_severe_consequence", "intergenic_variant"),
                impact="MODIFIER",
                found=True,
            )
        return EnsemblVEPAnnotation(
            consequence_type=record.get("most_severe_consequence"),
            found=True,
        )

    # Take the first (most severe with pick=1) consequence
    csq = consequences[0]

    # Extract SIFT prediction and score
    sift_prediction = csq.get("sift_prediction")
    sift_score = csq.get("sift_score")

    # Extract PolyPhen prediction and score
    polyphen_prediction = csq.get("polyphen_prediction")
    polyphen_score = csq.get("polyphen_score")

    # Extract protein domain
    domains = csq.get("domains", [])
    protein_domain = None
    if domains:
        # Prefer Pfam or InterPro domains
        for domain in domains:
            if domain.get("db", "").lower() in ("pfam", "interpro"):
                protein_domain = domain.get("name", domain.get("db"))
                break
        if not protein_domain and domains:
            protein_domain = domains[0].get("name", domains[0].get("db"))

    return EnsemblVEPAnnotation(
        consequence_type=csq.get("consequence_terms", [None])[0] if csq.get("consequence_terms") else record.get("most_severe_consequence"),
        impact=csq.get("impact"),
        gene_symbol=csq.get("gene_symbol"),
        gene_id=csq.get("gene_id"),
        transcript_id=csq.get("transcript_id"),
        biotype=csq.get("biotype"),
        amino_acid_change=csq.get("amino_acids"),
        codon_change=csq.get("codons"),
        sift_prediction=sift_prediction,
        sift_score=float(sift_score) if sift_score is not None else None,
        polyphen_prediction=polyphen_prediction,
        polyphen_score=float(polyphen_score) if polyphen_score is not None else None,
        protein_domain=protein_domain,
        exon=csq.get("exon"),
        found=True,
    )


async def query_vep(
    variant: Variant,
    client: httpx.AsyncClient | None = None,
) -> EnsemblVEPAnnotation:
    """Query Ensembl VEP for variant consequence prediction.

    Args:
        variant: The variant to annotate.
        client: Optional httpx.AsyncClient for connection reuse.

    Returns:
        EnsemblVEPAnnotation with consequence, impact, SIFT/PolyPhen scores.
        Returns found=False if the query fails.
    """
    url = _build_vep_url(variant)
    params = {
        "content-type": "application/json",
        "pick": "1",
        "SIFT": "b",
        "PolyPhen": "b",
        "domains": "1",
        "canonical": "1",
        "hgvs": "1",
    }
    headers = {"Content-Type": "application/json"}

    logger.info("Querying Ensembl VEP: %s", url)

    should_close = client is None
    if client is None:
        client = httpx.AsyncClient(timeout=30.0)

    try:
        response = await client.get(url, params=params, headers=headers)

        # Handle rate limiting
        if response.status_code == 429:
            retry_after = response.headers.get("Retry-After", "1")
            logger.warning("VEP rate limited, retry after %s seconds", retry_after)
            return EnsemblVEPAnnotation(found=False)

        response.raise_for_status()
        data = response.json()

        annotation = _parse_vep_response(data)
        logger.info(
            "VEP result for %s: %s (%s)",
            variant.variant_id,
            annotation.consequence_type,
            annotation.impact,
        )
        return annotation

    except httpx.HTTPStatusError as e:
        logger.error("VEP HTTP error: %s", e)
        return EnsemblVEPAnnotation(found=False)
    except Exception as e:
        logger.error("VEP query failed: %s", e)
        return EnsemblVEPAnnotation(found=False)
    finally:
        if should_close:
            await client.aclose()


async def query_vep_safe(
    variant: Variant,
    client: httpx.AsyncClient | None = None,
) -> tuple[EnsemblVEPAnnotation, str | None]:
    """Query VEP with error capture — never raises.

    Returns:
        Tuple of (annotation, error_message). error_message is None on success.
    """
    try:
        annotation = await query_vep(variant, client)
        return annotation, None
    except Exception as e:
        return EnsemblVEPAnnotation(found=False), f"VEP error: {e}"
