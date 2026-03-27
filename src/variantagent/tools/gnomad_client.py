"""gnomAD GraphQL API client.

Queries gnomAD for population allele frequencies via their GraphQL endpoint.

Critical notes:
- gnomAD has NO REST API — GraphQL only (POST to /api)
- The `af` field does NOT exist on VariantPopulation — must compute from ac/an
- Undocumented rate limit of ~10 requests/minute (6-second delay recommended)
- The API will silently block your IP without warning headers
"""

from __future__ import annotations

import asyncio
import logging

import httpx

from variantagent.models.annotation import GnomADFrequency
from variantagent.models.variant import Variant

logger = logging.getLogger(__name__)

GNOMAD_API = "https://gnomad.broadinstitute.org/api"

# Conservative delay between requests to avoid IP blocking
GNOMAD_DELAY_SECONDS = 6.0

VARIANT_QUERY = """
query GnomadVariant($variantId: String!, $datasetId: DatasetId!) {
  variant(variantId: $variantId, dataset: $datasetId) {
    exome {
      ac
      an
      ac_hom
      populations {
        id
        ac
        an
      }
    }
    genome {
      ac
      an
      ac_hom
      populations {
        id
        ac
        an
      }
    }
  }
}
"""


def _build_variant_id(variant: Variant) -> str:
    """Build gnomAD-style variant ID: {chrom}-{pos}-{ref}-{alt} (no 'chr' prefix)."""
    chrom = variant.normalized_chromosome
    return f"{chrom}-{variant.position}-{variant.reference}-{variant.alternate}"


def _compute_af(ac: int | None, an: int | None) -> float | None:
    """Compute allele frequency from allele count and allele number."""
    if ac is None or an is None or an == 0:
        return None
    return ac / an


def _parse_gnomad_response(data: dict) -> GnomADFrequency:
    """Parse gnomAD GraphQL response into GnomADFrequency.

    Merges exome and genome data when both are available.
    """
    variant_data = data.get("data", {}).get("variant")
    if variant_data is None:
        return GnomADFrequency(found=False)

    exome = variant_data.get("exome")
    genome = variant_data.get("genome")

    if exome is None and genome is None:
        return GnomADFrequency(found=False)

    # Use exome data preferentially, fall back to genome
    primary = exome or genome

    total_ac = (exome.get("ac", 0) if exome else 0) + (genome.get("ac", 0) if genome else 0)
    total_an = (exome.get("an", 0) if exome else 0) + (genome.get("an", 0) if genome else 0)
    total_hom = (exome.get("ac_hom", 0) if exome else 0) + (genome.get("ac_hom", 0) if genome else 0)

    # Build population frequency map from primary dataset
    pop_map: dict[str, dict[str, int]] = {}
    for source in [exome, genome]:
        if source is None:
            continue
        for pop in source.get("populations", []):
            pop_id = pop.get("id", "").lower()
            if pop_id not in pop_map:
                pop_map[pop_id] = {"ac": 0, "an": 0}
            pop_map[pop_id]["ac"] += pop.get("ac", 0)
            pop_map[pop_id]["an"] += pop.get("an", 0)

    def pop_af(pop_id: str) -> float | None:
        if pop_id in pop_map:
            return _compute_af(pop_map[pop_id]["ac"], pop_map[pop_id]["an"])
        return None

    return GnomADFrequency(
        overall_af=_compute_af(total_ac, total_an),
        afr_af=pop_af("afr"),
        amr_af=pop_af("amr"),
        asj_af=pop_af("asj"),
        eas_af=pop_af("eas"),
        fin_af=pop_af("fin"),
        nfe_af=pop_af("nfe"),
        sas_af=pop_af("sas"),
        homozygote_count=total_hom,
        allele_count=total_ac,
        allele_number=total_an,
        found=True,
    )


async def query_gnomad(
    variant: Variant,
    client: httpx.AsyncClient | None = None,
    dataset: str = "gnomad_r4",
) -> GnomADFrequency:
    """Query gnomAD for population allele frequencies.

    Args:
        variant: The variant to look up.
        client: Optional httpx.AsyncClient for connection reuse.
        dataset: gnomAD dataset version (default: gnomad_r4).

    Returns:
        GnomADFrequency with population-specific allele frequencies.
        Returns found=False if variant is not in gnomAD.
    """
    variant_id = _build_variant_id(variant)
    logger.info("Querying gnomAD: %s (dataset: %s)", variant_id, dataset)

    should_close = client is None
    if client is None:
        client = httpx.AsyncClient(timeout=30.0)

    try:
        payload = {
            "query": VARIANT_QUERY,
            "variables": {
                "variantId": variant_id,
                "datasetId": dataset,
            },
        }

        response = await client.post(
            GNOMAD_API,
            json=payload,
            headers={"Content-Type": "application/json"},
        )
        response.raise_for_status()
        data = response.json()

        # Check for GraphQL errors
        if "errors" in data:
            logger.warning("gnomAD GraphQL errors: %s", data["errors"])
            return GnomADFrequency(found=False)

        result = _parse_gnomad_response(data)

        if result.found:
            logger.info(
                "gnomAD result for %s: AF=%s",
                variant_id,
                f"{result.overall_af:.6f}" if result.overall_af is not None else "N/A",
            )
        else:
            logger.info("Variant not found in gnomAD: %s", variant_id)

        # Respect gnomAD's aggressive rate limiting
        await asyncio.sleep(GNOMAD_DELAY_SECONDS)

        return result

    except httpx.HTTPStatusError as e:
        logger.error("gnomAD HTTP error: %s", e)
        return GnomADFrequency(found=False)
    except Exception as e:
        logger.error("gnomAD query failed: %s", e)
        return GnomADFrequency(found=False)
    finally:
        if should_close:
            await client.aclose()


async def query_gnomad_safe(
    variant: Variant,
    client: httpx.AsyncClient | None = None,
    dataset: str = "gnomad_r4",
) -> tuple[GnomADFrequency, str | None]:
    """Query gnomAD with error capture — never raises.

    Returns:
        Tuple of (frequency_data, error_message). error_message is None on success.
    """
    try:
        result = await query_gnomad(variant, client, dataset)
        return result, None
    except Exception as e:
        return GnomADFrequency(found=False), f"gnomAD error: {e}"
