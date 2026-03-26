"""Reviewer Agent: Cross-checks conclusions and detects potential hallucinations.

This agent critiques the outputs of all other agents. It has NO external data
access — it works purely on the internal evidence already gathered.

System Prompt Role: Skeptical reviewer who challenges the classification.
Tools: claim_extractor, source_verifier, contradiction_detector
Distinct Because: Only agent that critiques other agents' outputs. No external API access.
"""

# TODO: Implement reviewer agent with:
# 1. Extract factual claims from the classification reasoning
# 2. Verify each claim against the annotation data (not the LLM's memory)
# 3. Detect contradictions between QC, annotation, and classification
#    e.g., "QC Agent says coverage is adequate but Classification Agent assumed low coverage"
# 4. Flag unsupported claims with hallucination risk level
# 5. Return ReviewerFindings for the TriageReport
