"""Annotation Agent: Queries public databases for variant annotation.

This agent is the sole interface to external APIs. It handles rate limiting,
retries, fallbacks, and error recovery for all database queries.

System Prompt Role: Database annotation specialist.
Tools (via MCP): clinvar_query, gnomad_frequency, ensembl_vep, uniprot_protein_impact
Distinct Because: Only agent interacting with external APIs. Handles networking concerns.
"""

# TODO: Implement annotation agent with:
# 1. ClinVar query via MCP server (or direct NCBI E-utilities)
# 2. gnomAD query via GraphQL API
# 3. Ensembl VEP via REST API
# 4. UniProt protein impact via REST API
# 5. Rate limiting and retry logic with exponential backoff
# 6. Graceful degradation: if one API fails, continue with others
# 7. Return VariantAnnotation with annotation_errors for any failures
