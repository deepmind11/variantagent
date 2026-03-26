"""Classification Agent: Applies ACMG/AMP criteria with chain-of-thought reasoning.

This agent combines evidence from QC, annotation, and literature to apply
ACMG/AMP variant classification criteria. It uses a deterministic rule engine
for the combining logic (which criteria combination = which classification)
while the LLM reasons about whether each individual criterion is met.

System Prompt Role: Variant classification specialist following ACMG/AMP 2015 standards.
Tools: acmg_rule_engine, confidence_scorer
Distinct Because: The only agent with a deterministic rule engine. LLM reasons about criteria;
    rule engine enforces correct combining logic.
"""

# TODO: Implement classification agent with:
# 1. Evaluate each ACMG criterion against available evidence
# 2. Chain-of-thought reasoning for each criterion decision
# 3. Deterministic ACMG combining rules (see acmg_engine.py in tools/)
# 4. Confidence scoring based on evidence strength and completeness
# 5. Return ACMGClassification with full criteria assessment
