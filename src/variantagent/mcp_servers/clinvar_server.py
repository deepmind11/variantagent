"""ClinVar MCP Server — reusable MCP server for querying NCBI ClinVar.

This server wraps the NCBI E-utilities API to provide ClinVar variant
lookups via the Model Context Protocol (MCP).

Can be used standalone or as part of the VariantAgent system.

Usage as standalone MCP server:
    python -m variantagent.mcp_servers.clinvar_server
"""

# TODO: Implement MCP server with:
# 1. Tool: clinvar_lookup
#    - Input: variant in HGVS notation, rsID, or genomic coordinates
#    - Output: ClinVarAnnotation (clinical significance, review status, conditions)
#    - Uses NCBI E-utilities (esearch + efetch)
#
# 2. Tool: clinvar_batch_lookup
#    - Input: list of variants
#    - Output: list of ClinVarAnnotations
#    - Batches API calls for efficiency
#
# 3. Rate limiting: respect NCBI's 3 req/sec (10/sec with API key)
# 4. Caching: LRU cache for repeated lookups
# 5. Error handling: graceful degradation if NCBI is down
