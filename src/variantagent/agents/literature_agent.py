"""Literature Agent: Searches published literature and ACMG guidelines via RAG.

This agent searches PubMed for variant-specific evidence and queries a
vector store containing embedded ACMG/AMP guidelines for relevant criteria.

System Prompt Role: Scientific literature specialist.
Tools: pubmed_search, rag_acmg_guidelines, citation_formatter
Distinct Because: Only agent using RAG. Combines live search + embedded knowledge base.
"""

# TODO: Implement literature agent with:
# 1. PubMed search via NCBI E-utilities (targeted queries for gene + variant + disease)
# 2. RAG over embedded ACMG/AMP 2015 guidelines (ChromaDB + sentence-transformers)
# 3. Citation formatting (PMID → formatted reference)
# 4. Relevance filtering: only return genuinely relevant publications
# 5. Summary generation: synthesize findings across multiple papers
