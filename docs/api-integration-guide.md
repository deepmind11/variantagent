# VariantAgent API Integration Guide

Practical implementation reference for the four public biological database APIs
that the annotation agent and literature agent depend on.

Last updated: 2026-03-27

---

## Table of Contents

1. [ClinVar (NCBI E-utilities)](#1-clinvar-ncbi-e-utilities)
2. [gnomAD (GraphQL API)](#2-gnomad-graphql-api)
3. [Ensembl VEP (REST API)](#3-ensembl-vep-rest-api)
4. [PubMed (NCBI E-utilities)](#4-pubmed-ncbi-e-utilities)
5. [Shared Infrastructure](#5-shared-infrastructure)
6. [References](#6-references)

---

## 1. ClinVar (NCBI E-utilities)

### Overview

ClinVar is NCBI's public archive of human genetic variants and their clinical
significance. Access is through the Entrez E-utilities HTTP API, using the
two-step `esearch` -> `esummary` (or `efetch`) pattern.

### Endpoints

| Endpoint   | URL                                                              | Purpose                        |
|------------|------------------------------------------------------------------|--------------------------------|
| esearch    | `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi`    | Search for variant IDs         |
| esummary   | `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi`   | Get structured variant summary |
| efetch     | `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi`     | Get full XML records           |

### Authentication

- **No API key**: 3 requests/second limit.
- **With API key**: 10 requests/second limit.
- Get a key from: https://www.ncbi.nlm.nih.gov/account/settings/ (NCBI account required).
- Pass as query parameter: `&api_key=YOUR_KEY`
- Always include `tool=variantagent` and `email=YOUR_EMAIL` parameters.

### Rate Limits

| Condition       | Limit             | Enforcement          |
|-----------------|-------------------|----------------------|
| No API key      | 3 req/sec per IP  | IP block if exceeded |
| With API key    | 10 req/sec per IP | IP block if exceeded |

NCBI does not return rate-limit headers. You must self-throttle.

### Query Formats

ClinVar esearch accepts the same query syntax as the web interface:

```
# By HGVS notation
NM_000546.6:c.743G>A

# By rsID
rs28934576[rsid]

# By genomic coordinates (gnomAD-style)
13:32932018:G:A

# By gene + significance
BRCA1[gene] AND pathogenic[clinsig]

# By variation ID
65533[uid]
```

### Response Structure (esummary JSON)

The esummary endpoint with `retmode=json` returns the most useful structured
data. Key fields in `result[uid]`:

```json
{
  "result": {
    "uids": ["65533"],
    "65533": {
      "uid": "65533",
      "accession": "VCV000065533.11",
      "title": "NM_000081.4(LYST):c.1540C>T (p.Arg514Ter)",
      "variation_set": [...],
      "germline_classification": {
        "description": "Pathogenic",
        "last_evaluated": "2025/10/18",
        "review_status": "criteria provided, multiple submitters, no conflicts",
        "trait_set": [
          {
            "trait_name": "Chediak-Higashi syndrome",
            "trait_xrefs": [
              {"db_source": "OMIM", "db_id": "214500"},
              {"db_source": "MedGen", "db_id": "C0007965"},
              {"db_source": "Orphanet", "db_id": "167"}
            ]
          }
        ]
      },
      "variation_type": "single nucleotide variant",
      "cdna_change": "c.1540C>T",
      "protein_change": "R514*",
      "molecular_consequence": "nonsense",
      "genes": [
        {"symbol": "LYST", "geneid": 1130, "strand": "-"}
      ],
      "chr": "1",
      "assembly_set": [
        {"assembly_name": "GRCh38", "chr": "1", "start": 235809278},
        {"assembly_name": "GRCh37", "chr": "1", "start": 235972578}
      ],
      "supporting_submissions": {
        "scv": ["SCV000092045", "SCV000778491", "SCV001245622", "SCV004874700"],
        "rcv": ["RCV000058510"]
      }
    }
  }
}
```

### Review Status to Stars Mapping

```python
REVIEW_STATUS_STARS: dict[str, int] = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, single submitter": 1,
    "no assertion criteria provided": 0,
    "no assertion provided": 0,
    "no classification provided": 0,
}
```

### Python Implementation

```python
"""ClinVar client using NCBI E-utilities."""

from __future__ import annotations

import asyncio
import xml.etree.ElementTree as ET
from typing import Any

import httpx

from variantagent.config import settings
from variantagent.models.annotation import ClinVarAnnotation

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

REVIEW_STATUS_STARS: dict[str, int] = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, single submitter": 1,
    "no assertion criteria provided": 0,
    "no assertion provided": 0,
    "no classification provided": 0,
}


def _base_params() -> dict[str, str]:
    """Return params that must accompany every E-utilities request."""
    params: dict[str, str] = {
        "tool": "variantagent",
        "email": settings.ncbi_email,
    }
    if settings.ncbi_api_key:
        params["api_key"] = settings.ncbi_api_key
    return params


async def _throttled_get(
    client: httpx.AsyncClient,
    url: str,
    params: dict[str, str],
    semaphore: asyncio.Semaphore,
) -> httpx.Response:
    """GET with concurrency-limited throttling."""
    async with semaphore:
        response = await client.get(url, params=params, timeout=30.0)
        response.raise_for_status()
        # Respect NCBI rate limit: 10/sec with key, 3/sec without
        delay = 0.1 if settings.ncbi_api_key else 0.34
        await asyncio.sleep(delay)
        return response


async def search_clinvar(
    client: httpx.AsyncClient,
    query: str,
    semaphore: asyncio.Semaphore,
    retmax: int = 5,
) -> list[str]:
    """Search ClinVar and return matching variation UIDs.

    Args:
        client: Reusable httpx async client.
        query: Search term (HGVS, rsID, coordinates, free text).
        semaphore: Concurrency limiter for rate control.
        retmax: Maximum results to return.

    Returns:
        List of ClinVar variation UIDs (strings).
    """
    params = {
        **_base_params(),
        "db": "clinvar",
        "term": query,
        "retmode": "json",
        "retmax": str(retmax),
    }
    response = await _throttled_get(
        client, f"{EUTILS_BASE}/esearch.fcgi", params, semaphore
    )
    data = response.json()
    return data.get("esearchresult", {}).get("idlist", [])


async def fetch_clinvar_summary(
    client: httpx.AsyncClient,
    uids: list[str],
    semaphore: asyncio.Semaphore,
) -> dict[str, Any]:
    """Fetch esummary JSON for a list of ClinVar UIDs.

    Args:
        client: Reusable httpx async client.
        uids: List of ClinVar variation UIDs from esearch.
        semaphore: Concurrency limiter for rate control.

    Returns:
        Raw esummary result dict keyed by UID.
    """
    if not uids:
        return {}
    params = {
        **_base_params(),
        "db": "clinvar",
        "id": ",".join(uids),
        "retmode": "json",
    }
    response = await _throttled_get(
        client, f"{EUTILS_BASE}/esummary.fcgi", params, semaphore
    )
    data = response.json()
    return data.get("result", {})


def _parse_clinvar_summary(uid: str, result: dict[str, Any]) -> ClinVarAnnotation:
    """Parse a single ClinVar esummary record into our model."""
    record = result.get(uid, {})
    if not record:
        return ClinVarAnnotation(found=False)

    # Extract germline classification (most common for constitutional variants)
    germline = record.get("germline_classification", {})
    clinical_sig = germline.get("description", None)
    review_status = germline.get("review_status", None)

    # Extract conditions from trait_set
    conditions: list[str] = []
    for trait_set in germline.get("trait_set", []):
        name = trait_set.get("trait_name")
        if name:
            conditions.append(name)

    # Count submitters from supporting_submissions
    scv_list = record.get("supporting_submissions", {}).get("scv", [])
    submitter_count = len(scv_list) if scv_list else None

    return ClinVarAnnotation(
        variation_id=record.get("accession"),
        clinical_significance=clinical_sig,
        review_status=review_status,
        review_stars=REVIEW_STATUS_STARS.get(review_status or "", None),
        conditions=conditions,
        submitter_count=submitter_count,
        last_evaluated=germline.get("last_evaluated"),
        found=True,
    )


async def query_clinvar(
    client: httpx.AsyncClient,
    query: str,
    semaphore: asyncio.Semaphore,
) -> ClinVarAnnotation:
    """Full pipeline: search ClinVar then fetch and parse the top result.

    Args:
        client: Reusable httpx async client.
        query: Variant identifier (HGVS, rsID, or coordinates).
        semaphore: Concurrency limiter for rate control.

    Returns:
        Parsed ClinVarAnnotation. If not found, found=False.

    Raises:
        httpx.HTTPStatusError: On non-retryable HTTP errors.
    """
    uids = await search_clinvar(client, query, semaphore)
    if not uids:
        return ClinVarAnnotation(found=False)

    result = await fetch_clinvar_summary(client, uids, semaphore)
    # Parse the first (most relevant) result
    return _parse_clinvar_summary(uids[0], result)
```

### Query Strategies by Input Type

```python
def build_clinvar_query(
    *,
    hgvs: str | None = None,
    rsid: str | None = None,
    chrom: str | None = None,
    pos: int | None = None,
    ref: str | None = None,
    alt: str | None = None,
) -> str:
    """Build the best ClinVar search term from available identifiers.

    Priority: HGVS > rsID > genomic coordinates.
    """
    if hgvs:
        # HGVS works directly as a search term
        return hgvs
    if rsid:
        # Bracket the rsID for field-specific search
        return f"{rsid}[rsid]" if not rsid.endswith("[rsid]") else rsid
    if chrom and pos and ref and alt:
        # gnomAD-style coordinate format
        return f"{chrom}:{pos}:{ref}:{alt}"
    msg = "At least one of hgvs, rsid, or (chrom, pos, ref, alt) is required"
    raise ValueError(msg)
```

### Error Handling

```python
async def query_clinvar_safe(
    client: httpx.AsyncClient,
    query: str,
    semaphore: asyncio.Semaphore,
    max_retries: int = 3,
) -> tuple[ClinVarAnnotation, str | None]:
    """Query ClinVar with retry logic and error capture.

    Returns:
        Tuple of (annotation, error_message). error_message is None on success.
    """
    last_error: str | None = None
    for attempt in range(max_retries):
        try:
            result = await query_clinvar(client, query, semaphore)
            return result, None
        except httpx.HTTPStatusError as exc:
            if exc.response.status_code == 429:
                # Rate limited -- back off exponentially
                wait = 2 ** attempt
                await asyncio.sleep(wait)
                last_error = f"ClinVar rate limited (attempt {attempt + 1})"
                continue
            if exc.response.status_code >= 500:
                wait = 2 ** attempt
                await asyncio.sleep(wait)
                last_error = f"ClinVar server error {exc.response.status_code}"
                continue
            # Client error (4xx) -- do not retry
            last_error = f"ClinVar HTTP {exc.response.status_code}: {exc.response.text[:200]}"
            break
        except httpx.TimeoutException:
            wait = 2 ** attempt
            await asyncio.sleep(wait)
            last_error = f"ClinVar timeout (attempt {attempt + 1})"
            continue
        except Exception as exc:
            last_error = f"ClinVar unexpected error: {exc!s}"
            break

    return ClinVarAnnotation(found=False), last_error
```

---

## 2. gnomAD (GraphQL API)

### Overview

gnomAD (Genome Aggregation Database) provides population allele frequencies
from 807,162 individuals (v4.1). The API is **GraphQL only** -- there is no
REST endpoint for variant queries. The API is designed to support the gnomAD
browser frontend; batch queries of arbitrary variant lists are not officially
supported.

### Endpoint

| Item        | Value                                        |
|-------------|----------------------------------------------|
| URL         | `https://gnomad.broadinstitute.org/api`      |
| Method      | POST                                         |
| Content-Type| `application/json`                           |
| Auth        | None required                                |

### Rate Limits

gnomAD has **undocumented rate limits** enforced at the IP level:

| Limit                           | Value              |
|---------------------------------|--------------------|
| Approximate safe rate           | 10 requests/minute |
| Recommended delay between calls | 6 seconds          |
| Blocking threshold              | ~10 rapid requests |

The API will silently block your IP if you exceed the rate limit. There are no
rate-limit response headers. The gnomAD team has stated the API is optimized
for browser use, not batch processing.

### Variant ID Format

gnomAD uses a dash-separated format: `{chrom}-{pos}-{ref}-{alt}`

```
# SNV on chromosome 17
17-7675088-C-T

# Insertion
2-233760233-C-CAT

# Deletion
1-55516888-GAAAC-G
```

### Available Datasets

| Dataset ID                | Description              |
|---------------------------|--------------------------|
| `gnomad_r4`               | gnomAD v4.1 (latest)    |
| `gnomad_r3`               | gnomAD v3.1.2           |
| `gnomad_r2_1`             | gnomAD v2.1.1           |
| `gnomad_r2_1_non_cancer`  | v2.1.1 non-cancer subset|
| `exac`                    | ExAC                    |

### Population IDs

| ID    | Population                   |
|-------|------------------------------|
| `afr` | African / African American   |
| `amr` | Latino / Admixed American    |
| `asj` | Ashkenazi Jewish             |
| `eas` | East Asian                   |
| `fin` | Finnish                      |
| `mid` | Middle Eastern (v4 only)     |
| `nfe` | Non-Finnish European         |
| `sas` | South Asian                  |
| `remaining` | Remaining / Other      |

### GraphQL Query Structure

```graphql
query GnomadVariant($variantId: String!, $datasetId: DatasetId!) {
  variant(variantId: $variantId, dataset: $datasetId) {
    variantId
    chrom
    pos
    ref
    alt
    flags
    exome {
      ac
      an
      ac_hom
      ac_hemi
      faf95 {
        popmax
        popmax_population
      }
      filters
      populations {
        id
        ac
        an
        ac_hom
        ac_hemi
      }
    }
    genome {
      ac
      an
      ac_hom
      ac_hemi
      faf95 {
        popmax
        popmax_population
      }
      filters
      populations {
        id
        ac
        an
        ac_hom
        ac_hemi
      }
    }
  }
}
```

**Important**: The `af` (allele frequency) field is NOT available on the
`populations` type. You must compute it: `af = ac / an` (when `an > 0`).

### Response Structure

```json
{
  "data": {
    "variant": {
      "variantId": "17-7675088-C-T",
      "chrom": "17",
      "pos": 7675088,
      "ref": "C",
      "alt": "T",
      "flags": [],
      "exome": {
        "ac": 5,
        "an": 1613150,
        "ac_hom": 0,
        "ac_hemi": 0,
        "faf95": {
          "popmax": 0.000008,
          "popmax_population": "nfe"
        },
        "filters": ["PASS"],
        "populations": [
          {"id": "afr", "ac": 0, "an": 123456, "ac_hom": 0, "ac_hemi": 0},
          {"id": "amr", "ac": 1, "an": 98765, "ac_hom": 0, "ac_hemi": 0},
          {"id": "asj", "ac": 0, "an": 25432, "ac_hom": 0, "ac_hemi": 0},
          {"id": "eas", "ac": 0, "an": 45678, "ac_hom": 0, "ac_hemi": 0},
          {"id": "fin", "ac": 0, "an": 67890, "ac_hom": 0, "ac_hemi": 0},
          {"id": "nfe", "ac": 4, "an": 890123, "ac_hom": 0, "ac_hemi": 0},
          {"id": "sas", "ac": 0, "an": 56789, "ac_hom": 0, "ac_hemi": 0}
        ]
      },
      "genome": null
    }
  }
}
```

Note: `genome` or `exome` can be `null` if data is unavailable for that
sequencing type. Always handle both being `null`.

### Python Implementation

```python
"""gnomAD client using the GraphQL API."""

from __future__ import annotations

import asyncio
from typing import Any

import httpx

from variantagent.models.annotation import GnomADFrequency

GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api"

# gnomAD rate limit: ~10 requests per minute. Use 6-second delays.
GNOMAD_REQUEST_DELAY = 6.0

GNOMAD_VARIANT_QUERY = """
query GnomadVariant($variantId: String!, $datasetId: DatasetId!) {
  variant(variantId: $variantId, dataset: $datasetId) {
    variantId
    chrom
    pos
    ref
    alt
    flags
    exome {
      ac
      an
      ac_hom
      filters
      populations {
        id
        ac
        an
        ac_hom
      }
    }
    genome {
      ac
      an
      ac_hom
      filters
      populations {
        id
        ac
        an
        ac_hom
      }
    }
  }
}
"""

# Mapping from gnomAD population ID to our model field name
_POP_FIELD_MAP: dict[str, str] = {
    "afr": "afr_af",
    "amr": "amr_af",
    "asj": "asj_af",
    "eas": "eas_af",
    "fin": "fin_af",
    "nfe": "nfe_af",
    "sas": "sas_af",
}


def _safe_af(ac: int | None, an: int | None) -> float | None:
    """Compute allele frequency, returning None if data is missing."""
    if ac is None or an is None or an == 0:
        return None
    return ac / an


def _merge_exome_genome(
    exome: dict[str, Any] | None,
    genome: dict[str, Any] | None,
) -> dict[str, Any]:
    """Merge exome and genome data, preferring combined counts.

    When both are available, sum ac and an for the overall frequency.
    For populations, sum per-population counts.
    """
    if exome is None and genome is None:
        return {}

    sources = [s for s in (exome, genome) if s is not None]

    total_ac = sum(s.get("ac", 0) or 0 for s in sources)
    total_an = sum(s.get("an", 0) or 0 for s in sources)
    total_hom = sum(s.get("ac_hom", 0) or 0 for s in sources)

    # Merge population data
    pop_data: dict[str, dict[str, int]] = {}
    for source in sources:
        for pop in source.get("populations", []):
            pid = pop["id"]
            if pid not in pop_data:
                pop_data[pid] = {"ac": 0, "an": 0, "ac_hom": 0}
            pop_data[pid]["ac"] += pop.get("ac", 0) or 0
            pop_data[pid]["an"] += pop.get("an", 0) or 0
            pop_data[pid]["ac_hom"] += pop.get("ac_hom", 0) or 0

    # Collect filter statuses
    filters: list[str] = []
    for source in sources:
        filters.extend(source.get("filters", []))

    return {
        "ac": total_ac,
        "an": total_an,
        "ac_hom": total_hom,
        "populations": pop_data,
        "filters": filters,
    }


def _parse_gnomad_response(data: dict[str, Any]) -> GnomADFrequency:
    """Parse gnomAD GraphQL response into our model."""
    variant = data.get("data", {}).get("variant")
    if variant is None:
        return GnomADFrequency(found=False)

    merged = _merge_exome_genome(variant.get("exome"), variant.get("genome"))
    if not merged:
        return GnomADFrequency(found=False)

    # Build population allele frequencies
    pop_afs: dict[str, float | None] = {}
    for pop_id, field_name in _POP_FIELD_MAP.items():
        pop = merged["populations"].get(pop_id, {})
        pop_afs[field_name] = _safe_af(pop.get("ac"), pop.get("an"))

    # Determine filtering status
    filters = merged.get("filters", [])
    filter_status = "PASS" if not filters or filters == ["PASS"] else ",".join(filters)

    return GnomADFrequency(
        overall_af=_safe_af(merged["ac"], merged["an"]),
        afr_af=pop_afs.get("afr_af"),
        amr_af=pop_afs.get("amr_af"),
        asj_af=pop_afs.get("asj_af"),
        eas_af=pop_afs.get("eas_af"),
        fin_af=pop_afs.get("fin_af"),
        nfe_af=pop_afs.get("nfe_af"),
        sas_af=pop_afs.get("sas_af"),
        homozygote_count=merged.get("ac_hom"),
        allele_count=merged.get("ac"),
        allele_number=merged.get("an"),
        filtering_status=filter_status,
        found=True,
    )


def build_variant_id(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Build gnomAD variant ID from components.

    Args:
        chrom: Chromosome (e.g., "17", "X"). No "chr" prefix.
        pos: 1-based genomic position.
        ref: Reference allele.
        alt: Alternate allele.

    Returns:
        gnomAD-format variant ID (e.g., "17-7675088-C-T").
    """
    # Strip "chr" prefix if present
    c = chrom.replace("chr", "")
    return f"{c}-{pos}-{ref}-{alt}"


async def query_gnomad(
    client: httpx.AsyncClient,
    variant_id: str,
    dataset: str = "gnomad_r4",
) -> GnomADFrequency:
    """Query gnomAD for a single variant's population frequencies.

    Args:
        client: Reusable httpx async client.
        variant_id: gnomAD-format ID (e.g., "17-7675088-C-T").
        dataset: Dataset identifier (default: gnomad_r4 for v4.1).

    Returns:
        Parsed GnomADFrequency. If not found, found=False.
    """
    payload = {
        "query": GNOMAD_VARIANT_QUERY,
        "variables": {
            "variantId": variant_id,
            "datasetId": dataset,
        },
    }
    response = await client.post(
        GNOMAD_API_URL,
        json=payload,
        headers={"Content-Type": "application/json"},
        timeout=30.0,
    )
    response.raise_for_status()
    data = response.json()

    # Check for GraphQL errors
    if "errors" in data:
        error_msgs = [e.get("message", "") for e in data["errors"]]
        raise ValueError(f"gnomAD GraphQL errors: {'; '.join(error_msgs)}")

    return _parse_gnomad_response(data)


async def query_gnomad_safe(
    client: httpx.AsyncClient,
    variant_id: str,
    dataset: str = "gnomad_r4",
    max_retries: int = 3,
) -> tuple[GnomADFrequency, str | None]:
    """Query gnomAD with retry logic and mandatory rate limiting.

    The 6-second delay between requests is ALWAYS applied to respect
    gnomAD's undocumented rate limits (~10 req/min).

    Returns:
        Tuple of (frequency_data, error_message).
    """
    last_error: str | None = None
    for attempt in range(max_retries):
        try:
            result = await query_gnomad(client, variant_id, dataset)
            # ALWAYS wait after a successful request
            await asyncio.sleep(GNOMAD_REQUEST_DELAY)
            return result, None
        except httpx.HTTPStatusError as exc:
            wait = GNOMAD_REQUEST_DELAY * (2 ** attempt)
            await asyncio.sleep(wait)
            last_error = f"gnomAD HTTP {exc.response.status_code}"
            continue
        except ValueError as exc:
            last_error = f"gnomAD query error: {exc!s}"
            break
        except httpx.TimeoutException:
            wait = GNOMAD_REQUEST_DELAY * (2 ** attempt)
            await asyncio.sleep(wait)
            last_error = f"gnomAD timeout (attempt {attempt + 1})"
            continue
        except Exception as exc:
            last_error = f"gnomAD unexpected error: {exc!s}"
            break

    return GnomADFrequency(found=False), last_error
```

### Fallback Strategy

If the gnomAD GraphQL API is slow or blocked, consider these alternatives:

1. **gnomAD Hail Table downloads**: Download precomputed VCF/TSV files from
   `https://gnomad.broadinstitute.org/downloads` and query locally.
2. **gnomAD SQLite (gnomad-db)**: PyPI package `gnomad-db` provides SQLite-based
   local queries for GRCh37/38 (requires ~30GB disk).
3. **Ensembl VEP colocated_variants**: The VEP response includes gnomAD frequencies
   in the `colocated_variants` array, avoiding a separate gnomAD call entirely.

---

## 3. Ensembl VEP (REST API)

### Overview

The Ensembl Variant Effect Predictor (VEP) determines the effect of variants
on genes, transcripts, and protein sequence. The REST API supports individual
GET queries and batch POST queries.

### Endpoints

| Endpoint          | Method | URL Pattern                                                   | Use Case                    |
|-------------------|--------|---------------------------------------------------------------|-----------------------------|
| VEP by HGVS      | GET    | `/vep/{species}/hgvs/{hgvs_notation}`                        | Single variant, HGVS input  |
| VEP by HGVS      | POST   | `/vep/{species}/hgvs`                                        | Batch (up to 200 variants)  |
| VEP by region     | GET    | `/vep/{species}/region/{region}/{allele}`                     | Single variant, coordinates |
| VEP by region     | POST   | `/vep/{species}/region`                                      | Batch by coordinates        |
| VEP by ID         | GET    | `/vep/{species}/id/{id}`                                     | Query by rsID               |

Base URL: `https://rest.ensembl.org`

### Rate Limits

| Metric            | Value                |
|-------------------|----------------------|
| Requests/hour     | 55,000               |
| Avg requests/sec  | ~15                  |
| POST batch size   | 200 variants max     |

Rate limit state is communicated via response headers:

```
X-RateLimit-Limit: 55000
X-RateLimit-Period: 3600
X-RateLimit-Remaining: 54999
X-RateLimit-Reset: 3599
```

When exceeded, the API returns HTTP 429 with a `Retry-After` header
(floating-point seconds).

### Query Formats

**GET by HGVS notation:**
```
GET /vep/human/hgvs/ENST00000366667:c.803C>T?content-type=application/json
GET /vep/human/hgvs/NM_000546.6:c.743G>A?content-type=application/json
GET /vep/human/hgvs/9:g.22125504G>C?content-type=application/json
```

**GET by genomic region:**
```
GET /vep/human/region/17:7675088-7675088:1/T?content-type=application/json
```
Region format: `{chrom}:{start}-{end}:{strand}/{allele}`

**POST batch by HGVS:**
```json
POST /vep/human/hgvs
Content-Type: application/json

{
  "hgvs_notations": [
    "ENST00000366667:c.803C>T",
    "NM_000546.6:c.743G>A"
  ]
}
```

**POST batch by region:**
```json
POST /vep/human/region
Content-Type: application/json

{
  "variants": [
    "17 7675088 7675088 C/T 1",
    "9 22125503 22125504 G/C 1"
  ]
}
```

### Recommended Query Parameters

For VariantAgent, always include these parameters to get full annotation:

```python
VEP_DEFAULT_PARAMS: dict[str, str] = {
    "content-type": "application/json",
    "SIFT": "b",           # Include SIFT prediction + score
    "PolyPhen": "b",       # Include PolyPhen prediction + score
    "CADD": "1",           # CADD deleteriousness score
    "canonical": "1",      # Flag canonical transcript
    "domains": "1",        # Protein domains
    "hgvs": "1",           # HGVS nomenclature
    "numbers": "1",        # Exon/intron numbers
    "mane": "1",           # MANE Select transcript
    "pick": "1",           # One consequence per variant (simplifies parsing)
    "variant_class": "1",  # Variant class (SNV, insertion, etc.)
}
```

### Response Structure

```json
[
  {
    "input": "ENST00000366667:c.803C>T",
    "assembly_name": "GRCh38",
    "seq_region_name": "1",
    "start": 230710021,
    "end": 230710021,
    "strand": 1,
    "allele_string": "C/T",
    "most_severe_consequence": "missense_variant",
    "transcript_consequences": [
      {
        "transcript_id": "ENST00000366667",
        "gene_id": "ENSG00000135744",
        "gene_symbol": "AGT",
        "biotype": "protein_coding",
        "canonical": 1,
        "consequence_terms": ["missense_variant"],
        "impact": "MODERATE",
        "hgvsc": "ENST00000366667.5:c.803C>T",
        "hgvsp": "ENSP00000355627.5:p.Ala268Val",
        "amino_acids": "A/V",
        "codons": "gCc/gTc",
        "protein_start": 268,
        "protein_end": 268,
        "exon": "3/5",
        "sift_prediction": "tolerated",
        "sift_score": 0.08,
        "polyphen_prediction": "benign",
        "polyphen_score": 0.003,
        "domains": [
          {"db": "Pfam", "name": "PF00079"},
          {"db": "PANTHER", "name": "PTHR11588"},
          {"db": "Gene3D", "name": "1.10.1370.10"}
        ]
      }
    ],
    "colocated_variants": [
      {
        "id": "rs1228544607",
        "allele_string": "C/T",
        "frequencies": {
          "T": {
            "gnomade": 0.00000899,
            "gnomade_nfe": 0.00000899
          }
        }
      }
    ]
  }
]
```

### Python Implementation

```python
"""Ensembl VEP client using the REST API."""

from __future__ import annotations

import asyncio
from typing import Any

import httpx

from variantagent.models.annotation import EnsemblVEPAnnotation

ENSEMBL_BASE_URL = "https://rest.ensembl.org"

VEP_DEFAULT_PARAMS: dict[str, str] = {
    "content-type": "application/json",
    "SIFT": "b",
    "PolyPhen": "b",
    "canonical": "1",
    "domains": "1",
    "hgvs": "1",
    "numbers": "1",
    "mane": "1",
    "pick": "1",
    "variant_class": "1",
}


async def _ensembl_get(
    client: httpx.AsyncClient,
    path: str,
    params: dict[str, str] | None = None,
) -> Any:
    """Make a rate-limit-aware GET request to the Ensembl REST API.

    Reads X-RateLimit-Remaining and Retry-After headers and adjusts
    pacing accordingly.
    """
    merged_params = {**VEP_DEFAULT_PARAMS, **(params or {})}
    response = await client.get(
        f"{ENSEMBL_BASE_URL}{path}",
        params=merged_params,
        timeout=30.0,
    )

    # Handle rate limiting
    if response.status_code == 429:
        retry_after = float(response.headers.get("Retry-After", "1.0"))
        await asyncio.sleep(retry_after)
        # Retry once after waiting
        response = await client.get(
            f"{ENSEMBL_BASE_URL}{path}",
            params=merged_params,
            timeout=30.0,
        )

    response.raise_for_status()

    # Adaptive throttle: slow down when nearing limit
    remaining = int(response.headers.get("X-RateLimit-Remaining", "55000"))
    if remaining < 1000:
        await asyncio.sleep(1.0)
    elif remaining < 5000:
        await asyncio.sleep(0.2)

    return response.json()


async def _ensembl_post(
    client: httpx.AsyncClient,
    path: str,
    body: dict[str, Any],
) -> Any:
    """Make a rate-limit-aware POST request to the Ensembl REST API."""
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    response = await client.post(
        f"{ENSEMBL_BASE_URL}{path}",
        json=body,
        headers=headers,
        timeout=60.0,  # POST can be slower
    )

    if response.status_code == 429:
        retry_after = float(response.headers.get("Retry-After", "2.0"))
        await asyncio.sleep(retry_after)
        response = await client.post(
            f"{ENSEMBL_BASE_URL}{path}",
            json=body,
            headers=headers,
            timeout=60.0,
        )

    response.raise_for_status()
    return response.json()


def _pick_canonical_consequence(data: dict[str, Any]) -> dict[str, Any] | None:
    """Select the canonical transcript consequence from VEP output.

    Priority: canonical transcript > MANE Select > first protein_coding.
    If pick=1 was used, there will typically be only one consequence.
    """
    consequences = data.get("transcript_consequences", [])
    if not consequences:
        return None

    # If only one consequence (pick=1 mode), return it
    if len(consequences) == 1:
        return consequences[0]

    # Prefer canonical
    for tc in consequences:
        if tc.get("canonical") == 1 and tc.get("biotype") == "protein_coding":
            return tc

    # Fallback: first protein_coding
    for tc in consequences:
        if tc.get("biotype") == "protein_coding":
            return tc

    return consequences[0]


def _parse_vep_response(data: list[dict[str, Any]]) -> EnsemblVEPAnnotation:
    """Parse VEP response (always returns a list) into our model."""
    if not data:
        return EnsemblVEPAnnotation(found=False)

    record = data[0]  # First variant result
    tc = _pick_canonical_consequence(record)

    if tc is None:
        return EnsemblVEPAnnotation(
            consequence_type=record.get("most_severe_consequence"),
            found=True,
        )

    # Extract protein domain as human-readable string
    domains = tc.get("domains", [])
    domain_str = None
    if domains:
        # Prefer Pfam, then PANTHER, then first available
        for db_name in ("Pfam", "PANTHER", "Gene3D"):
            for d in domains:
                if d.get("db") == db_name:
                    domain_str = f"{d['db']}:{d['name']}"
                    break
            if domain_str:
                break
        if not domain_str:
            d = domains[0]
            domain_str = f"{d.get('db', 'unknown')}:{d.get('name', 'unknown')}"

    # Extract amino acid change in compact format (e.g., "R175H")
    aa_change = None
    amino_acids = tc.get("amino_acids", "")
    protein_start = tc.get("protein_start")
    if "/" in amino_acids and protein_start:
        ref_aa, alt_aa = amino_acids.split("/", 1)
        aa_change = f"{ref_aa}{protein_start}{alt_aa}"

    return EnsemblVEPAnnotation(
        consequence_type=record.get("most_severe_consequence"),
        impact=tc.get("impact"),
        gene_symbol=tc.get("gene_symbol"),
        gene_id=tc.get("gene_id"),
        transcript_id=tc.get("transcript_id"),
        biotype=tc.get("biotype"),
        amino_acid_change=aa_change,
        codon_change=tc.get("codons"),
        sift_prediction=tc.get("sift_prediction"),
        sift_score=tc.get("sift_score"),
        polyphen_prediction=tc.get("polyphen_prediction"),
        polyphen_score=tc.get("polyphen_score"),
        protein_domain=domain_str,
        exon=tc.get("exon"),
        found=True,
    )


async def query_vep_hgvs(
    client: httpx.AsyncClient,
    hgvs_notation: str,
) -> EnsemblVEPAnnotation:
    """Query VEP for a single variant using HGVS notation.

    Args:
        client: Reusable httpx async client.
        hgvs_notation: e.g., "NM_000546.6:c.743G>A" or "17:g.7675088C>T"

    Returns:
        Parsed VEP annotation.
    """
    # URL-encode the HGVS notation (the > character in particular)
    import urllib.parse
    encoded = urllib.parse.quote(hgvs_notation, safe="")
    path = f"/vep/human/hgvs/{encoded}"
    data = await _ensembl_get(client, path)
    return _parse_vep_response(data)


async def query_vep_region(
    client: httpx.AsyncClient,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
) -> EnsemblVEPAnnotation:
    """Query VEP for a single variant using genomic coordinates.

    Args:
        client: Reusable httpx async client.
        chrom: Chromosome (e.g., "17").
        pos: 1-based position.
        ref: Reference allele.
        alt: Alternate allele.
    """
    c = chrom.replace("chr", "")

    if len(ref) == 1 and len(alt) == 1:
        # SNV: region is pos-pos
        region = f"{c}:{pos}-{pos}:1"
        allele = alt
    elif len(ref) > len(alt):
        # Deletion
        start = pos + 1
        end = pos + len(ref) - 1
        region = f"{c}:{start}-{end}:1"
        allele = "-"
    else:
        # Insertion
        region = f"{c}:{pos}-{pos}:1"
        allele = alt[len(ref):]  # Inserted bases only

    path = f"/vep/human/region/{region}/{allele}"
    data = await _ensembl_get(client, path)
    return _parse_vep_response(data)


async def query_vep_batch(
    client: httpx.AsyncClient,
    hgvs_notations: list[str],
) -> list[EnsemblVEPAnnotation]:
    """Query VEP for a batch of variants (max 200 per request).

    Args:
        client: Reusable httpx async client.
        hgvs_notations: List of HGVS notation strings.

    Returns:
        List of parsed VEP annotations (one per input variant).
    """
    results: list[EnsemblVEPAnnotation] = []
    # Chunk into batches of 200
    for i in range(0, len(hgvs_notations), 200):
        chunk = hgvs_notations[i : i + 200]
        body = {"hgvs_notations": chunk}
        data = await _ensembl_post(client, "/vep/human/hgvs", body)
        # POST returns a flat list of results, one per input
        for record in data:
            annotation = _parse_vep_response([record])
            results.append(annotation)
    return results


async def query_vep_safe(
    client: httpx.AsyncClient,
    *,
    hgvs: str | None = None,
    chrom: str | None = None,
    pos: int | None = None,
    ref: str | None = None,
    alt: str | None = None,
    max_retries: int = 3,
) -> tuple[EnsemblVEPAnnotation, str | None]:
    """Query VEP with retry logic and error capture.

    Accepts either HGVS notation or genomic coordinates.

    Returns:
        Tuple of (annotation, error_message).
    """
    last_error: str | None = None
    for attempt in range(max_retries):
        try:
            if hgvs:
                result = await query_vep_hgvs(client, hgvs)
            elif chrom and pos and ref and alt:
                result = await query_vep_region(client, chrom, pos, ref, alt)
            else:
                return (
                    EnsemblVEPAnnotation(found=False),
                    "VEP requires either hgvs or (chrom, pos, ref, alt)",
                )
            return result, None
        except httpx.HTTPStatusError as exc:
            if exc.response.status_code == 429:
                retry_after = float(
                    exc.response.headers.get("Retry-After", str(2 ** attempt))
                )
                await asyncio.sleep(retry_after)
                last_error = f"VEP rate limited (attempt {attempt + 1})"
                continue
            if exc.response.status_code >= 500:
                await asyncio.sleep(2 ** attempt)
                last_error = f"VEP server error {exc.response.status_code}"
                continue
            if exc.response.status_code == 400:
                # Bad request -- likely invalid HGVS notation. Do not retry.
                last_error = f"VEP bad request: {exc.response.text[:200]}"
                break
            last_error = f"VEP HTTP {exc.response.status_code}"
            break
        except httpx.TimeoutException:
            await asyncio.sleep(2 ** attempt)
            last_error = f"VEP timeout (attempt {attempt + 1})"
            continue
        except Exception as exc:
            last_error = f"VEP unexpected error: {exc!s}"
            break

    return EnsemblVEPAnnotation(found=False), last_error
```

---

## 4. PubMed (NCBI E-utilities)

### Overview

PubMed literature search uses the same NCBI E-utilities infrastructure as
ClinVar. The workflow is `esearch` (get PMIDs) -> `efetch` (get abstracts)
or `esummary` (get citation metadata).

### Endpoints

Same base as ClinVar: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils`

| Endpoint   | Use                                        |
|------------|--------------------------------------------|
| esearch    | Search PubMed, get list of PMIDs           |
| efetch     | Get full abstracts (XML or text)           |
| esummary   | Get citation metadata (JSON supported)     |

### Authentication and Rate Limits

Same as ClinVar (shared NCBI infrastructure):

- **No API key**: 3 req/sec
- **With API key**: 10 req/sec
- Use the same `tool`, `email`, and `api_key` parameters.

### Search Strategy for Variants

Building effective PubMed queries for variant evidence:

```python
def build_pubmed_query(
    gene: str,
    variant: str | None = None,
    disease: str | None = None,
    hgvs: str | None = None,
) -> str:
    """Build a targeted PubMed query for variant literature evidence.

    Strategy: Start specific, broaden if no results.
    """
    parts: list[str] = []

    # Gene is always required
    parts.append(f"({gene}[Gene Name] OR {gene}[Title/Abstract])")

    # Variant identifiers -- try multiple representations
    if hgvs:
        parts.append(f'("{hgvs}"[Title/Abstract])')
    elif variant:
        # variant could be "R175H", "p.Arg175His", "rs28934576", etc.
        parts.append(f'("{variant}"[Title/Abstract])')

    # Disease context narrows results
    if disease:
        parts.append(f"({disease}[Title/Abstract])")

    return " AND ".join(parts)


# Example queries for TP53 R175H:
# Specific:  (TP53[Gene Name]) AND ("p.Arg175His"[Title/Abstract])
# Broader:   (TP53[Gene Name]) AND ("R175H"[Title/Abstract])
# With disease: (TP53[Gene Name]) AND ("R175H"[Title/Abstract]) AND (cancer[Title/Abstract])
```

### esearch Response (JSON)

```json
{
  "esearchresult": {
    "count": "42",
    "retmax": "20",
    "retstart": "0",
    "idlist": ["39876543", "38765432", "37654321", "..."],
    "querytranslation": "(TP53[Gene Name]) AND (R175H[Title/Abstract])"
  }
}
```

### efetch for Abstracts (XML)

Use `rettype=abstract&retmode=xml` to get structured PubMed XML:

```xml
<PubmedArticleSet>
  <PubmedArticle>
    <MedlineCitation>
      <PMID>39876543</PMID>
      <Article>
        <ArticleTitle>TP53 R175H gain-of-function mutations...</ArticleTitle>
        <Abstract>
          <AbstractText>Background: The TP53 p.Arg175His...</AbstractText>
        </Abstract>
        <AuthorList>
          <Author><LastName>Smith</LastName><Initials>JD</Initials></Author>
        </AuthorList>
        <Journal>
          <Title>Nature Genetics</Title>
          <JournalIssue>
            <Volume>55</Volume>
            <PubDate><Year>2024</Year></PubDate>
          </JournalIssue>
        </Journal>
      </Article>
    </MedlineCitation>
  </PubmedArticle>
</PubmedArticleSet>
```

### esummary for Citation Metadata (JSON)

```json
{
  "result": {
    "uids": ["39876543"],
    "39876543": {
      "uid": "39876543",
      "pubdate": "2024 Mar",
      "source": "Nat Genet",
      "authors": [
        {"name": "Smith JD", "authtype": "Author"},
        {"name": "Doe AB", "authtype": "Author"}
      ],
      "title": "TP53 R175H gain-of-function mutations...",
      "volume": "55",
      "issue": "3",
      "pages": "234-245",
      "fulljournalname": "Nature genetics",
      "elocationid": "doi: 10.1038/s41588-024-01234-5",
      "pubtype": ["Journal Article"]
    }
  }
}
```

### Python Implementation

```python
"""PubMed client using NCBI E-utilities."""

from __future__ import annotations

import asyncio
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from typing import Any

import httpx

from variantagent.config import settings

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


@dataclass(frozen=True)
class PubMedArticle:
    """Minimal PubMed article representation."""

    pmid: str
    title: str
    abstract: str
    authors: list[str]
    journal: str
    year: str
    doi: str | None = None

    def formatted_citation(self) -> str:
        """Return a formatted citation string."""
        author_str = ", ".join(self.authors[:3])
        if len(self.authors) > 3:
            author_str += ", et al."
        doi_str = f" doi:{self.doi}" if self.doi else ""
        return f"{author_str}. {self.title} {self.journal}. {self.year}.{doi_str} PMID:{self.pmid}"


def _base_params() -> dict[str, str]:
    """Return params that must accompany every E-utilities request."""
    params: dict[str, str] = {
        "tool": "variantagent",
        "email": settings.ncbi_email,
    }
    if settings.ncbi_api_key:
        params["api_key"] = settings.ncbi_api_key
    return params


async def search_pubmed(
    client: httpx.AsyncClient,
    query: str,
    semaphore: asyncio.Semaphore,
    retmax: int = 10,
    sort: str = "relevance",
) -> list[str]:
    """Search PubMed and return matching PMIDs.

    Args:
        client: Reusable httpx async client.
        query: PubMed search query.
        semaphore: Concurrency limiter for rate control.
        retmax: Maximum results (default 10, max 10000).
        sort: Sort order -- "relevance" (default), "pub_date", "Author".

    Returns:
        List of PMID strings.
    """
    params = {
        **_base_params(),
        "db": "pubmed",
        "term": query,
        "retmode": "json",
        "retmax": str(retmax),
        "sort": sort,
    }
    async with semaphore:
        response = await client.get(
            f"{EUTILS_BASE}/esearch.fcgi",
            params=params,
            timeout=30.0,
        )
        response.raise_for_status()
        delay = 0.1 if settings.ncbi_api_key else 0.34
        await asyncio.sleep(delay)

    data = response.json()
    return data.get("esearchresult", {}).get("idlist", [])


async def fetch_pubmed_summaries(
    client: httpx.AsyncClient,
    pmids: list[str],
    semaphore: asyncio.Semaphore,
) -> dict[str, Any]:
    """Fetch esummary JSON for a list of PMIDs.

    Returns:
        Raw esummary result dict keyed by PMID.
    """
    if not pmids:
        return {}
    params = {
        **_base_params(),
        "db": "pubmed",
        "id": ",".join(pmids),
        "retmode": "json",
    }
    async with semaphore:
        response = await client.get(
            f"{EUTILS_BASE}/esummary.fcgi",
            params=params,
            timeout=30.0,
        )
        response.raise_for_status()
        delay = 0.1 if settings.ncbi_api_key else 0.34
        await asyncio.sleep(delay)

    return response.json().get("result", {})


async def fetch_pubmed_abstracts(
    client: httpx.AsyncClient,
    pmids: list[str],
    semaphore: asyncio.Semaphore,
) -> list[PubMedArticle]:
    """Fetch full abstracts via efetch XML and parse into PubMedArticle objects.

    Args:
        client: Reusable httpx async client.
        pmids: List of PMID strings.
        semaphore: Concurrency limiter.

    Returns:
        List of PubMedArticle with title, abstract, authors, journal, year.
    """
    if not pmids:
        return []

    params = {
        **_base_params(),
        "db": "pubmed",
        "id": ",".join(pmids),
        "rettype": "abstract",
        "retmode": "xml",
    }
    async with semaphore:
        response = await client.get(
            f"{EUTILS_BASE}/efetch.fcgi",
            params=params,
            timeout=30.0,
        )
        response.raise_for_status()
        delay = 0.1 if settings.ncbi_api_key else 0.34
        await asyncio.sleep(delay)

    return _parse_pubmed_xml(response.text)


def _parse_pubmed_xml(xml_text: str) -> list[PubMedArticle]:
    """Parse PubMed efetch XML into PubMedArticle objects."""
    articles: list[PubMedArticle] = []
    root = ET.fromstring(xml_text)

    for article_elem in root.findall(".//PubmedArticle"):
        citation = article_elem.find("MedlineCitation")
        if citation is None:
            continue

        # PMID
        pmid_elem = citation.find("PMID")
        pmid = pmid_elem.text if pmid_elem is not None else "unknown"

        # Article details
        article = citation.find("Article")
        if article is None:
            continue

        # Title
        title_elem = article.find("ArticleTitle")
        title = title_elem.text or "" if title_elem is not None else ""

        # Abstract -- may have multiple AbstractText elements
        abstract_parts: list[str] = []
        abstract_elem = article.find("Abstract")
        if abstract_elem is not None:
            for at in abstract_elem.findall("AbstractText"):
                label = at.get("Label", "")
                text = at.text or ""
                if label:
                    abstract_parts.append(f"{label}: {text}")
                else:
                    abstract_parts.append(text)
        abstract = " ".join(abstract_parts)

        # Authors
        authors: list[str] = []
        author_list = article.find("AuthorList")
        if author_list is not None:
            for author in author_list.findall("Author"):
                last = author.find("LastName")
                initials = author.find("Initials")
                if last is not None:
                    name = last.text or ""
                    if initials is not None and initials.text:
                        name += f" {initials.text}"
                    authors.append(name)

        # Journal
        journal_elem = article.find("Journal/Title")
        journal = journal_elem.text or "" if journal_elem is not None else ""

        # Year
        year = ""
        pub_date = article.find("Journal/JournalIssue/PubDate")
        if pub_date is not None:
            year_elem = pub_date.find("Year")
            if year_elem is not None:
                year = year_elem.text or ""

        # DOI
        doi: str | None = None
        for eid in article.findall("ELocationID"):
            if eid.get("EIdType") == "doi":
                doi = eid.text

        articles.append(
            PubMedArticle(
                pmid=pmid,
                title=title,
                abstract=abstract,
                authors=authors,
                journal=journal,
                year=year,
                doi=doi,
            )
        )

    return articles


def build_pubmed_query(
    gene: str,
    variant: str | None = None,
    disease: str | None = None,
    hgvs: str | None = None,
) -> str:
    """Build a targeted PubMed query for variant literature evidence.

    Uses field-qualified terms for precision.
    """
    parts: list[str] = []
    parts.append(f"({gene}[Gene Name] OR {gene}[Title/Abstract])")

    if hgvs:
        parts.append(f'("{hgvs}"[Title/Abstract])')
    elif variant:
        parts.append(f'("{variant}"[Title/Abstract])')

    if disease:
        parts.append(f"({disease}[Title/Abstract])")

    return " AND ".join(parts)


async def search_variant_literature(
    client: httpx.AsyncClient,
    gene: str,
    semaphore: asyncio.Semaphore,
    variant: str | None = None,
    disease: str | None = None,
    hgvs: str | None = None,
    max_results: int = 5,
) -> tuple[list[PubMedArticle], str | None]:
    """Search PubMed for variant-relevant literature with progressive broadening.

    Strategy:
    1. Search with most specific query (gene + HGVS + disease).
    2. If no results, broaden (gene + variant protein change).
    3. If still no results, broadest (gene + disease).

    Returns:
        Tuple of (articles, error_message).
    """
    queries_to_try: list[str] = []

    # Most specific first
    if hgvs:
        queries_to_try.append(build_pubmed_query(gene, hgvs=hgvs, disease=disease))
    if variant:
        queries_to_try.append(build_pubmed_query(gene, variant=variant, disease=disease))
    if variant and not disease:
        queries_to_try.append(build_pubmed_query(gene, variant=variant))
    # Broadest fallback: gene + disease
    if disease:
        queries_to_try.append(build_pubmed_query(gene, disease=disease))

    # Deduplicate while preserving order
    seen: set[str] = set()
    unique_queries: list[str] = []
    for q in queries_to_try:
        if q not in seen:
            seen.add(q)
            unique_queries.append(q)

    for query in unique_queries:
        try:
            pmids = await search_pubmed(client, query, semaphore, retmax=max_results)
            if pmids:
                articles = await fetch_pubmed_abstracts(client, pmids, semaphore)
                return articles, None
        except Exception as exc:
            return [], f"PubMed error: {exc!s}"

    return [], None  # No results found, but no error either
```

---

## 5. Shared Infrastructure

### Async HTTP Client Factory

All four API clients share a single `httpx.AsyncClient` for connection pooling:

```python
"""Shared HTTP client and rate-limiting infrastructure."""

from __future__ import annotations

import asyncio
from contextlib import asynccontextmanager
from collections.abc import AsyncGenerator

import httpx

from variantagent.config import settings


@asynccontextmanager
async def create_api_client() -> AsyncGenerator[httpx.AsyncClient, None]:
    """Create a shared async HTTP client with connection pooling.

    Usage:
        async with create_api_client() as client:
            result = await query_clinvar(client, "NM_000546.6:c.743G>A", semaphore)
    """
    async with httpx.AsyncClient(
        follow_redirects=True,
        limits=httpx.Limits(
            max_connections=10,
            max_keepalive_connections=5,
        ),
    ) as client:
        yield client


def create_ncbi_semaphore() -> asyncio.Semaphore:
    """Create a semaphore for NCBI rate limiting.

    With API key: 10 concurrent requests (self-throttled to 10/sec).
    Without API key: 3 concurrent requests (self-throttled to 3/sec).
    """
    max_concurrent = 10 if settings.ncbi_api_key else 3
    return asyncio.Semaphore(max_concurrent)
```

### Orchestrating All Queries in the Annotation Agent

```python
"""Example: how the annotation agent orchestrates all API calls."""

from __future__ import annotations

import asyncio

import httpx

from variantagent.models.annotation import VariantAnnotation

# Import all clients (from code above)
# from variantagent.clients.clinvar import query_clinvar_safe, build_clinvar_query
# from variantagent.clients.gnomad import query_gnomad_safe, build_variant_id
# from variantagent.clients.ensembl_vep import query_vep_safe
# from variantagent.clients.pubmed import search_variant_literature
# from variantagent.clients.shared import create_api_client, create_ncbi_semaphore


async def annotate_variant(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    gene: str | None = None,
    hgvs: str | None = None,
) -> VariantAnnotation:
    """Annotate a single variant from all four databases.

    Queries ClinVar, Ensembl VEP, and PubMed in parallel.
    gnomAD runs sequentially due to its strict rate limits.
    Each query has independent error handling -- a failure in one
    database does not block the others.
    """
    async with create_api_client() as client:
        ncbi_semaphore = create_ncbi_semaphore()
        annotation_errors: list[str] = []

        # --- Parallel group: ClinVar + VEP + PubMed ---
        clinvar_query = build_clinvar_query(
            hgvs=hgvs, chrom=chrom, pos=pos, ref=ref, alt=alt
        )

        clinvar_task = query_clinvar_safe(client, clinvar_query, ncbi_semaphore)
        vep_task = query_vep_safe(
            client, hgvs=hgvs, chrom=chrom, pos=pos, ref=ref, alt=alt
        )

        # PubMed needs gene symbol; may come from VEP if not provided
        pubmed_articles: list[str] = []

        clinvar_result, vep_result = await asyncio.gather(
            clinvar_task, vep_task
        )

        clinvar_annotation, clinvar_error = clinvar_result
        vep_annotation, vep_error = vep_result

        if clinvar_error:
            annotation_errors.append(clinvar_error)
        if vep_error:
            annotation_errors.append(vep_error)

        # Use gene from VEP if not provided
        resolved_gene = gene or vep_annotation.gene_symbol

        # PubMed search (depends on gene symbol)
        if resolved_gene:
            variant_str = vep_annotation.amino_acid_change
            articles, pubmed_error = await search_variant_literature(
                client,
                resolved_gene,
                ncbi_semaphore,
                variant=variant_str,
            )
            pubmed_articles = [a.pmid for a in articles]
            if pubmed_error:
                annotation_errors.append(pubmed_error)

        # --- Sequential: gnomAD (strict rate limits) ---
        gnomad_id = build_variant_id(chrom, pos, ref, alt)
        gnomad_annotation, gnomad_error = await query_gnomad_safe(
            client, gnomad_id
        )
        if gnomad_error:
            annotation_errors.append(gnomad_error)

        return VariantAnnotation(
            clinvar=clinvar_annotation,
            gnomad=gnomad_annotation,
            ensembl_vep=vep_annotation,
            pubmed_references=pubmed_articles,
            annotation_errors=annotation_errors,
        )
```

### Caching Strategy

For repeated lookups (common in batch VCF processing), add an LRU cache:

```python
from functools import lru_cache
from typing import Any

# For synchronous wrappers or pre-computed lookups:
@lru_cache(maxsize=1024)
def cached_variant_key(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Normalize variant into a cache key."""
    return f"{chrom}-{pos}-{ref}-{alt}"


# For async, use a dict-based cache:
class AsyncVariantCache:
    """Simple async-safe variant annotation cache."""

    def __init__(self, maxsize: int = 1024) -> None:
        self._cache: dict[str, VariantAnnotation] = {}
        self._maxsize = maxsize

    def get(self, key: str) -> VariantAnnotation | None:
        return self._cache.get(key)

    def put(self, key: str, value: VariantAnnotation) -> None:
        if len(self._cache) >= self._maxsize:
            # Evict oldest entry (FIFO)
            oldest = next(iter(self._cache))
            del self._cache[oldest]
        self._cache[key] = value
```

---

## 6. References

### NCBI E-utilities

- E-utilities documentation: https://www.ncbi.nlm.nih.gov/books/NBK25497/
- Parameters reference: https://www.ncbi.nlm.nih.gov/books/NBK25499/
- Sample applications: https://www.ncbi.nlm.nih.gov/books/NBK25498/
- ClinVar access guide: https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/
- ClinVar FTP primer: https://www.ncbi.nlm.nih.gov/clinvar/docs/ftp_primer/
- API key registration: https://www.ncbi.nlm.nih.gov/account/settings/

### gnomAD

- gnomAD browser: https://gnomad.broadinstitute.org
- GraphQL endpoint: https://gnomad.broadinstitute.org/api
- gnomAD browser source: https://github.com/broadinstitute/gnomad-browser
- gnomAD v4 batch tool: https://github.com/sfbizzari/gnomADv4-Batch-tool-pythonAPI
- gnomAD rate limit discussion: https://discuss.gnomad.broadinstitute.org/t/blocked-when-using-api-to-get-af/149
- DeepWiki GraphQL API reference: https://deepwiki.com/broadinstitute/gnomad-browser/4-graphql-api

### Ensembl VEP

- REST API home: https://rest.ensembl.org
- VEP HGVS GET: https://rest.ensembl.org/documentation/info/vep_hgvs_get
- VEP HGVS POST: https://rest.ensembl.org/documentation/info/vep_hgvs_post
- VEP region GET: https://rest.ensembl.org/documentation/info/vep_region_get
- Rate limits: https://github.com/Ensembl/ensembl-rest/wiki/Rate-Limits
- POST requests: https://github.com/Ensembl/ensembl-rest/wiki/POST-Requests

### PubMed

- PubMed API guide: https://library.cumc.columbia.edu/kb/getting-started-pubmed-api
- Abstract batch download: https://github.com/erilu/pubmed-abstract-compiler

### Open-Source Projects with Relevant Implementations

- BRCA Exchange gnomAD pipeline: https://github.com/BRCAChallenge/brca-exchange/tree/master/pipeline/gnomad
- gnomAD Python API (deprecated but instructive): https://github.com/furkanmtorun/gnomad_python_api
- gnomAD batch query gist: https://gist.github.com/ressy/6fd7f6ee6401ac8e703dc2709399869e
- metapub ClinVar fetcher: https://metapub.readthedocs.io/en/latest/_modules/metapub/clinvarfetcher.html
- pyEnsemblRest: https://github.com/gawbul/pyEnsemblRest
