# Phase 2 Implementation Plan

## Pre-Phase: Fix Code Review Issues (30 min)

Before adding features, fix the 4 HIGH issues from code review:

### Fix 1: Add LangGraph reducers to list fields in AnalysisState
- Add `Annotated[list[...], operator.add]` to `provenance`, `errors`, `reviewer_findings`
- Simplify all nodes to return `{"provenance": [new_entry]}` instead of `state["provenance"] + [new_entry]`
- This prevents silent data loss if we ever add parallel branches

### Fix 2: ACMG engine — add missing Likely Benign rule
- Add `2x Supporting Benign = Likely Benign` rule (Richards Table 5)
- Add exhaustive boundary tests for all combining rules

### Fix 3: QC fail report — say "not assessed" instead of "Unknown"
- When classification is None (QC abort), generate appropriate summary
- "Variant interpretation was not performed due to QC failure" instead of "classified as Unknown"

### Fix 4: QC agent — report all likely causes, not just first
- Change `QCIssue.likely_cause: str` → `likely_causes: list[str]`
- Update all `assess_flagstat` and `assess_multiqc` calls
- Update tests

**Commit:** `fix: resolve code review issues (reducers, ACMG rules, QC reporting)`

---

## Step 1: Annotation Agent — ClinVar Client (2-3 hours)

**What:** Implement the ClinVar query using NCBI E-utilities.

**Implementation:**
- `src/variantagent/tools/clinvar_client.py`
  - `async def query_clinvar(variant: Variant) -> ClinVarAnnotation`
  - Two-step: `esearch` (get UIDs) → `esummary` (get JSON data)
  - Query formats: rsID (`rs28934576[rsid]`), HGVS, genomic coords
  - Parse: `germline_classification.description`, `review_status`, `trait_set`
  - Rate limiting: `asyncio.Semaphore(3)` (or 10 with API key)
  - Retry with exponential backoff on 429/500
  - Return `ClinVarAnnotation` with `found=True/False`

**Test strategy:**
- Record real API responses with VCR.py for known variants (TP53 R175H, BRCA1 185delAG)
- Unit test parser logic against recorded responses
- One live integration test (marked `@pytest.mark.live`, skipped in CI)

**Commit:** `feat: implement ClinVar client with NCBI E-utilities`

---

## Step 2: Annotation Agent — Ensembl VEP Client (1-2 hours)

**What:** Query Ensembl VEP REST API for variant consequence prediction.

**Implementation:**
- `src/variantagent/tools/ensembl_client.py`
  - `async def query_vep(variant: Variant) -> EnsemblVEPAnnotation`
  - GET `/vep/human/region/{chr}:{pos}:{pos}/{alt}?pick=1&SIFT=b&PolyPhen=b&domains=1&canonical=1`
  - Parse: consequence_type, impact, SIFT/PolyPhen scores, protein domain, amino acid change
  - Rate limit: generous (15/sec), but respect `Retry-After` headers
  - Return `EnsemblVEPAnnotation` with `found=True/False`

**Test strategy:** Same VCR.py pattern as ClinVar.

**Commit:** `feat: implement Ensembl VEP client`

---

## Step 3: Annotation Agent — gnomAD Client (1-2 hours)

**What:** Query gnomAD GraphQL API for population allele frequencies.

**Implementation:**
- `src/variantagent/tools/gnomad_client.py`
  - `async def query_gnomad(variant: Variant) -> GnomADFrequency`
  - GraphQL POST to `https://gnomad.broadinstitute.org/api`
  - Variant ID format: `{chrom}-{pos}-{ref}-{alt}` (no "chr" prefix)
  - Compute AF from `ac/an` (the `af` field doesn't exist on `VariantPopulation`)
  - Handle: exome and/or genome fields can be null
  - Rate limit: 6-second delay between requests (aggressive gnomAD limit)
  - Return `GnomADFrequency` with population-specific AFs

**Risk:** gnomAD's undocumented rate limiting can silently block your IP. Use conservative delays.

**Test strategy:** VCR.py cassettes. Mock the 6s delay in tests.

**Commit:** `feat: implement gnomAD GraphQL client`

---

## Step 4: Wire Annotation Agent into Orchestrator (1 hour)

**What:** Replace the placeholder `annotation_node` with real API calls.

**Implementation:**
- Update `annotation_node` in `orchestrator.py` to call all three clients
- Run ClinVar + VEP in parallel (both are fast), gnomAD sequentially
- Collect errors from each — if one fails, continue with others
- Update provenance entries with actual data sources queried

**Test strategy:** Integration test with VCR.py cassettes for the full annotation flow.

**Commit:** `feat: wire annotation agent to live ClinVar, gnomAD, VEP APIs`

---

## Step 5: Literature Agent — PubMed Search (1-2 hours)

**What:** Search PubMed for variant-specific evidence.

**Implementation:**
- `src/variantagent/tools/pubmed_client.py`
  - `async def search_pubmed(gene: str, variant: str, max_results: int = 5) -> list[dict]`
  - Progressive broadening: gene + HGVS → gene + protein change → gene + disease
  - `esearch` → `esummary` (JSON) for citation metadata
  - Return PMIDs, titles, journals, years

- Update `literature_node` in orchestrator to call PubMed client

**Commit:** `feat: implement PubMed search client and literature agent`

---

## Step 6: Literature Agent — RAG over ACMG Guidelines (2-3 hours)

**What:** Embed ACMG/AMP 2015 guidelines into ChromaDB for semantic search.

**Implementation:**
- `src/variantagent/tools/rag.py`
  - `build_knowledge_base()` — chunk and embed ACMG guidelines PDF/text
  - `query_knowledge_base(question: str, k: int = 3) -> list[str]` — semantic search
  - Use `sentence-transformers/all-MiniLM-L6-v2` for embeddings
  - ChromaDB for storage (embedded, no server needed)

- `data/knowledge_base/acmg_guidelines.md` — ACMG criteria descriptions (manually curated from the open-access paper)

- Update `literature_node` to query both PubMed and ACMG knowledge base

**Commit:** `feat: add RAG knowledge base with ACMG guidelines`

---

## Step 7: Classification Agent — LLM Criterion Assessment (2-3 hours)

**What:** Use the LLM to evaluate which ACMG criteria are met, then feed to the deterministic rule engine.

**Implementation:**
- Update `classification_node` in orchestrator:
  - Build a prompt with: variant info, ClinVar annotation, gnomAD frequencies, VEP consequence, literature findings, ACMG criteria descriptions (from RAG)
  - Ask the LLM: "For each of these ACMG criteria, is it met? Provide reasoning."
  - Parse LLM output into `ACMGCriteria` (structured output / function calling)
  - Feed criteria to `acmg_engine.classify()` for deterministic classification
  - Return `ACMGClassification` with full reasoning chain

**Key design:** The LLM does the judgment. The rule engine does the math. This is the core architectural insight.

**Test strategy:**
- Mock the LLM response with a known set of criteria
- Verify the rule engine produces the correct classification
- Test with a real LLM call for known pathogenic variant (TP53 R175H) — marked `@pytest.mark.live`

**Commit:** `feat: implement LLM-based ACMG criterion assessment`

---

## Step 8: Reviewer Agent — Self-Evaluation (1-2 hours)

**What:** Expand the reviewer to do real claim extraction and verification.

**Implementation:**
- Expand `review_node` in orchestrator:
  - Extract factual claims from classification reasoning
  - Cross-check: does ClinVar agree with the classification?
  - Cross-check: does the population frequency (gnomAD) align with the pathogenicity call?
  - Flag contradictions (e.g., "classified as Pathogenic but gnomAD AF > 1%")
  - Assign hallucination risk scores to each claim
  - Compute overall confidence from evidence completeness + consistency

**Commit:** `feat: implement reviewer agent with claim verification`

---

## Step 9: CLI and Demo (1-2 hours)

**What:** Wire the CLI to run end-to-end with rich terminal output.

**Implementation:**
- Update `cli.py`:
  - Parse variant input (string format or VCF path)
  - Call `analyze_variant()`
  - Display results with `rich` (colored panels, tables, provenance tree)
  - JSON output option (`--output report.json`)

**Commit:** `feat: wire CLI for end-to-end variant analysis`

---

## Step 10: Polish and Push (1-2 hours)

- Update README with actual demo output
- Run full test suite, fix any failures
- Ensure CI passes (lint, typecheck, test)
- Final commit and push

**Commit:** `docs: update README with real demo output and usage examples`

---

## Implementation Order (Dependencies)

```
Pre-Phase (fixes) ─────────────────────────────────────────►
    │
    ├── Step 1: ClinVar client ──┐
    ├── Step 2: VEP client ──────┤── Step 4: Wire annotation agent
    ├── Step 3: gnomAD client ───┘         │
    │                                      │
    ├── Step 5: PubMed client ─────────────┤
    ├── Step 6: RAG knowledge base ────────┤
    │                                      │
    │                               Step 7: Classification agent
    │                                      │
    │                               Step 8: Reviewer agent
    │                                      │
    │                               Step 9: CLI
    │                                      │
    │                               Step 10: Polish
```

Steps 1-3 can be built in parallel (independent API clients).
Steps 5-6 can be built in parallel with steps 1-3.
Steps 7-10 are sequential (each depends on the previous).

---

## Git Commit Strategy (10 natural commits)

```
fix: resolve code review issues (reducers, ACMG rules, QC reporting)
feat: implement ClinVar client with NCBI E-utilities
feat: implement Ensembl VEP client
feat: implement gnomAD GraphQL client
feat: wire annotation agent to live ClinVar, gnomAD, VEP APIs
feat: implement PubMed search client and literature agent
feat: add RAG knowledge base with ACMG guidelines
feat: implement LLM-based ACMG criterion assessment
feat: implement reviewer agent with claim verification
feat: wire CLI for end-to-end variant analysis
```

---

## Risk Assessment

| Risk | Impact | Mitigation |
|------|--------|-----------|
| gnomAD rate limits block IP | Can't demo | Use VCR.py cassettes; cache aggressively; fallback to VEP's colocated_variants |
| LLM structured output parsing fails | Wrong classification | Use function calling / tool_use for structured output; validate with Pydantic |
| NCBI E-utilities downtime | Annotation incomplete | Graceful degradation — annotate with whatever succeeds |
| ChromaDB embedding quality poor for bio text | RAG returns irrelevant chunks | Use domain-appropriate chunking; test retrieval quality manually |
| Cost of LLM calls during development | Budget | Use Claude Haiku / GPT-4o-mini for dev; only test with Sonnet for benchmarks |

---

## Estimated Total: 15-20 hours

Realistically 2-3 focused sessions to complete all of Phase 2.
