"""Deterministic ACMG/AMP variant classification combining rules.

This implements the combining logic from Richards et al. 2015 (Table 5).
The LLM agents decide WHICH criteria are met; this engine determines
the final classification by applying the rule table.

This is NOT an LLM — it is a deterministic rule engine that enforces
correct classification logic. This separation ensures reproducibility.
"""

from variantagent.models.classification import (
    ACMGClassificationResult,
    ACMGCriteria,
    EvidenceDirection,
    EvidenceStrength,
)


def classify(criteria: ACMGCriteria) -> tuple[ACMGClassificationResult, str]:
    """Apply ACMG/AMP combining rules to determine classification.

    Implements Table 5 from Richards et al. 2015:
    https://doi.org/10.1038/gim.2015.30

    Args:
        criteria: Evaluated ACMG criteria with applied/not-applied flags.

    Returns:
        Tuple of (classification result, rule description that produced it).
    """
    pathogenic_codes = criteria.get_pathogenic_codes()
    benign_codes = criteria.get_benign_codes()

    # Count by strength
    pvs = sum(1 for c in pathogenic_codes if c.strength == EvidenceStrength.VERY_STRONG)
    ps = sum(1 for c in pathogenic_codes if c.strength == EvidenceStrength.STRONG)
    pm = sum(1 for c in pathogenic_codes if c.strength == EvidenceStrength.MODERATE)
    pp = sum(1 for c in pathogenic_codes if c.strength == EvidenceStrength.SUPPORTING)

    ba = sum(1 for c in benign_codes if c.strength == EvidenceStrength.VERY_STRONG)
    bs = sum(1 for c in benign_codes if c.strength == EvidenceStrength.STRONG)
    bp = sum(1 for c in benign_codes if c.strength == EvidenceStrength.SUPPORTING)

    # Check for conflicting evidence
    has_pathogenic = len(pathogenic_codes) > 0
    has_benign = len(benign_codes) > 0

    # === BENIGN RULES (check first — BA1 is standalone) ===

    # BA1 alone = Benign (standalone rule)
    if ba >= 1:
        return ACMGClassificationResult.BENIGN, "BA1 (standalone benign: allele frequency > 5%)"

    # >= 2 Strong Benign = Benign
    if bs >= 2:
        return ACMGClassificationResult.BENIGN, f"{bs} Strong Benign criteria"

    # Likely Benign: 1 Strong + 1 Supporting
    if bs >= 1 and bp >= 1:
        return ACMGClassificationResult.LIKELY_BENIGN, "1 Strong Benign + 1 Supporting Benign"

    # === PATHOGENIC RULES ===

    # Pathogenic (i): 1 Very Strong + >= 1 Strong
    if pvs >= 1 and ps >= 1:
        return ACMGClassificationResult.PATHOGENIC, f"1 Very Strong + {ps} Strong"

    # Pathogenic (ii): 1 Very Strong + >= 2 Moderate
    if pvs >= 1 and pm >= 2:
        return ACMGClassificationResult.PATHOGENIC, f"1 Very Strong + {pm} Moderate"

    # Pathogenic (iii): 1 Very Strong + 1 Moderate + 1 Supporting
    if pvs >= 1 and pm >= 1 and pp >= 1:
        return ACMGClassificationResult.PATHOGENIC, "1 Very Strong + 1 Moderate + 1 Supporting"

    # Pathogenic (iv): 1 Very Strong + >= 2 Supporting
    if pvs >= 1 and pp >= 2:
        return ACMGClassificationResult.PATHOGENIC, f"1 Very Strong + {pp} Supporting"

    # Pathogenic (v): >= 2 Strong
    if ps >= 2:
        return ACMGClassificationResult.PATHOGENIC, f"{ps} Strong"

    # Pathogenic (vi): 1 Strong + >= 3 Moderate
    if ps >= 1 and pm >= 3:
        return ACMGClassificationResult.PATHOGENIC, f"1 Strong + {pm} Moderate"

    # Pathogenic (vii): 1 Strong + 2 Moderate + >= 2 Supporting
    if ps >= 1 and pm >= 2 and pp >= 2:
        return ACMGClassificationResult.PATHOGENIC, f"1 Strong + {pm} Moderate + {pp} Supporting"

    # Pathogenic (viii): 1 Strong + 1 Moderate + >= 4 Supporting
    if ps >= 1 and pm >= 1 and pp >= 4:
        return ACMGClassificationResult.PATHOGENIC, f"1 Strong + 1 Moderate + {pp} Supporting"

    # === LIKELY PATHOGENIC RULES ===

    # Likely Pathogenic (i): 1 Very Strong + 1 Moderate
    if pvs >= 1 and pm == 1:
        return ACMGClassificationResult.LIKELY_PATHOGENIC, "1 Very Strong + 1 Moderate"

    # Likely Pathogenic (ii): 1 Strong + 1-2 Moderate
    if ps >= 1 and 1 <= pm <= 2:
        return ACMGClassificationResult.LIKELY_PATHOGENIC, f"1 Strong + {pm} Moderate"

    # Likely Pathogenic (iii): 1 Strong + >= 2 Supporting
    if ps >= 1 and pp >= 2:
        return ACMGClassificationResult.LIKELY_PATHOGENIC, f"1 Strong + {pp} Supporting"

    # Likely Pathogenic (iv): >= 3 Moderate
    if pm >= 3:
        return ACMGClassificationResult.LIKELY_PATHOGENIC, f"{pm} Moderate"

    # Likely Pathogenic (v): 2 Moderate + >= 2 Supporting
    if pm >= 2 and pp >= 2:
        return ACMGClassificationResult.LIKELY_PATHOGENIC, f"{pm} Moderate + {pp} Supporting"

    # Likely Pathogenic (vi): 1 Moderate + >= 4 Supporting
    if pm >= 1 and pp >= 4:
        return ACMGClassificationResult.LIKELY_PATHOGENIC, f"1 Moderate + {pp} Supporting"

    # === VUS: anything that doesn't meet the above rules ===
    if has_pathogenic and has_benign:
        return ACMGClassificationResult.VUS, "Conflicting pathogenic and benign evidence"

    if has_pathogenic:
        return ACMGClassificationResult.VUS, "Insufficient pathogenic evidence for classification"

    if has_benign:
        return ACMGClassificationResult.LIKELY_BENIGN, "Benign evidence present but insufficient for definitive classification"

    return ACMGClassificationResult.VUS, "No evidence criteria met"
