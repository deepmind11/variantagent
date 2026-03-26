# ADR-003: Deterministic ACMG Rule Engine Separate from LLM

## Status
Accepted

## Context
ACMG/AMP variant classification involves two distinct steps:
1. Deciding which evidence criteria are met (requires judgment)
2. Combining criteria into a classification (deterministic rules from Table 5 of Richards et al. 2015)

## Decision
Separate these into: **LLM-powered criterion assessment** (Classification Agent) + **deterministic rule engine** (`acmg_engine.py`).

## Rationale
1. **Reproducibility.** Given the same set of applied criteria, the classification will always be the same. The LLM may vary in which criteria it applies, but the combining logic is guaranteed correct.
2. **Auditability.** In regulated environments, the classification logic must be inspectable and verifiable. A deterministic rule engine can be unit-tested exhaustively.
3. **Trust.** Separating "judgment" from "rules" makes the system more trustworthy. Users can verify the rule engine independently of the LLM.
4. **Testing.** The rule engine has ~20 distinct combining rules that can be exhaustively tested. LLM outputs cannot.

## Consequences
- The Classification Agent must produce structured output (which criteria are met/not met) rather than a free-text classification.
- The rule engine must be kept in sync with ACMG/AMP standards if they are updated.
- Users can override individual criterion assessments and re-run the rule engine for a different classification.
