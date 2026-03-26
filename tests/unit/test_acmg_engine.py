"""Tests for the deterministic ACMG combining rules engine."""

from variantagent.models.classification import (
    ACMGClassificationResult,
    ACMGCriteria,
    EvidenceCode,
    EvidenceDirection,
    EvidenceStrength,
)
from variantagent.tools.acmg_engine import classify


def _make_code(
    code: str,
    direction: EvidenceDirection,
    strength: EvidenceStrength,
    applied: bool = True,
) -> EvidenceCode:
    """Helper to create an EvidenceCode for testing."""
    return EvidenceCode(
        code=code,
        name=f"Test {code}",
        direction=direction,
        strength=strength,
        applied=applied,
        reasoning="Test reasoning",
        data_source="test",
        confidence=0.9,
    )


class TestACMGEngine:
    def test_ba1_standalone_benign(self) -> None:
        """BA1 alone (AF > 5%) should classify as Benign."""
        criteria = ACMGCriteria(
            ba1=_make_code("BA1", EvidenceDirection.BENIGN, EvidenceStrength.VERY_STRONG)
        )
        result, rule = classify(criteria)
        assert result == ACMGClassificationResult.BENIGN
        assert "BA1" in rule

    def test_pathogenic_pvs1_plus_ps1(self) -> None:
        """1 Very Strong + 1 Strong = Pathogenic."""
        criteria = ACMGCriteria(
            pvs1=_make_code("PVS1", EvidenceDirection.PATHOGENIC, EvidenceStrength.VERY_STRONG),
            ps1=_make_code("PS1", EvidenceDirection.PATHOGENIC, EvidenceStrength.STRONG),
        )
        result, rule = classify(criteria)
        assert result == ACMGClassificationResult.PATHOGENIC

    def test_likely_pathogenic_pvs1_plus_pm2(self) -> None:
        """1 Very Strong + 1 Moderate = Likely Pathogenic."""
        criteria = ACMGCriteria(
            pvs1=_make_code("PVS1", EvidenceDirection.PATHOGENIC, EvidenceStrength.VERY_STRONG),
            pm2=_make_code("PM2", EvidenceDirection.PATHOGENIC, EvidenceStrength.MODERATE),
        )
        result, rule = classify(criteria)
        assert result == ACMGClassificationResult.LIKELY_PATHOGENIC

    def test_likely_benign_bs_plus_bp(self) -> None:
        """1 Strong Benign + 1 Supporting Benign = Likely Benign."""
        criteria = ACMGCriteria(
            bs1=_make_code("BS1", EvidenceDirection.BENIGN, EvidenceStrength.STRONG),
            bp4=_make_code("BP4", EvidenceDirection.BENIGN, EvidenceStrength.SUPPORTING),
        )
        result, rule = classify(criteria)
        assert result == ACMGClassificationResult.LIKELY_BENIGN

    def test_vus_no_evidence(self) -> None:
        """No evidence = VUS."""
        criteria = ACMGCriteria()
        result, rule = classify(criteria)
        assert result == ACMGClassificationResult.VUS

    def test_vus_conflicting_evidence(self) -> None:
        """Pathogenic + Benign evidence without meeting any threshold = VUS."""
        criteria = ACMGCriteria(
            pp3=_make_code("PP3", EvidenceDirection.PATHOGENIC, EvidenceStrength.SUPPORTING),
            bp4=_make_code("BP4", EvidenceDirection.BENIGN, EvidenceStrength.SUPPORTING),
        )
        result, rule = classify(criteria)
        assert result == ACMGClassificationResult.VUS

    def test_unapplied_codes_ignored(self) -> None:
        """Codes with applied=False should not count."""
        criteria = ACMGCriteria(
            pvs1=_make_code(
                "PVS1", EvidenceDirection.PATHOGENIC, EvidenceStrength.VERY_STRONG, applied=False
            ),
            ps1=_make_code(
                "PS1", EvidenceDirection.PATHOGENIC, EvidenceStrength.STRONG, applied=False
            ),
        )
        result, rule = classify(criteria)
        assert result == ACMGClassificationResult.VUS
