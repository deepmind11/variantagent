"""Tests for the classification agent's evidence-based criterion evaluation."""

from variantagent.agents.orchestrator import _evaluate_criteria_from_evidence, _calculate_confidence
from variantagent.models.annotation import (
    ClinVarAnnotation, EnsemblVEPAnnotation, GnomADFrequency, VariantAnnotation,
)
from variantagent.models.variant import Variant, VariantType
from variantagent.tools.acmg_engine import classify


class TestEvaluateCriteriaFromEvidence:
    def _make_variant(self) -> Variant:
        return Variant(
            chromosome="chr17", position=7674220, reference="G", alternate="A",
            gene="TP53", variant_type=VariantType.SNV,
        )

    def test_common_variant_gets_ba1(self) -> None:
        """AF > 5% → BA1 (standalone benign)."""
        annotation = VariantAnnotation(
            gnomad=GnomADFrequency(found=True, overall_af=0.15),
        )
        criteria = _evaluate_criteria_from_evidence(self._make_variant(), annotation)
        assert criteria.ba1 is not None
        assert criteria.ba1.applied is True
        assert criteria.ba1.code == "BA1"

        result, _ = classify(criteria)
        assert result.value == "Benign"

    def test_rare_variant_gets_pm2(self) -> None:
        """AF < 0.0001 → PM2 (absent from population databases)."""
        annotation = VariantAnnotation(
            gnomad=GnomADFrequency(found=True, overall_af=0.00001),
        )
        criteria = _evaluate_criteria_from_evidence(self._make_variant(), annotation)
        assert criteria.pm2 is not None
        assert criteria.pm2.applied is True

    def test_absent_variant_gets_pm2(self) -> None:
        """Not found in gnomAD → PM2."""
        annotation = VariantAnnotation(
            gnomad=GnomADFrequency(found=False),
        )
        criteria = _evaluate_criteria_from_evidence(self._make_variant(), annotation)
        assert criteria.pm2 is not None
        assert criteria.pm2.applied is True

    def test_clinvar_pathogenic_gets_pp5(self) -> None:
        """ClinVar Pathogenic with 2+ stars → PP5."""
        annotation = VariantAnnotation(
            clinvar=ClinVarAnnotation(
                found=True, clinical_significance="Pathogenic",
                review_stars=3, submitter_count=10,
            ),
        )
        criteria = _evaluate_criteria_from_evidence(self._make_variant(), annotation)
        assert criteria.pp5 is not None
        assert criteria.pp5.applied is True

    def test_clinvar_benign_gets_bp6(self) -> None:
        """ClinVar Benign with 2+ stars → BP6."""
        annotation = VariantAnnotation(
            clinvar=ClinVarAnnotation(
                found=True, clinical_significance="Benign",
                review_stars=2,
            ),
        )
        criteria = _evaluate_criteria_from_evidence(self._make_variant(), annotation)
        assert criteria.bp6 is not None
        assert criteria.bp6.applied is True

    def test_deleterious_computational_gets_pp3(self) -> None:
        """SIFT deleterious + PolyPhen damaging → PP3."""
        annotation = VariantAnnotation(
            ensembl_vep=EnsemblVEPAnnotation(
                found=True, sift_prediction="deleterious", sift_score=0.0,
                polyphen_prediction="probably_damaging", polyphen_score=0.999,
            ),
        )
        criteria = _evaluate_criteria_from_evidence(self._make_variant(), annotation)
        assert criteria.pp3 is not None
        assert criteria.pp3.applied is True

    def test_benign_computational_gets_bp4(self) -> None:
        """SIFT tolerated + PolyPhen benign → BP4."""
        annotation = VariantAnnotation(
            ensembl_vep=EnsemblVEPAnnotation(
                found=True, sift_prediction="tolerated", sift_score=0.8,
                polyphen_prediction="benign", polyphen_score=0.01,
            ),
        )
        criteria = _evaluate_criteria_from_evidence(self._make_variant(), annotation)
        assert criteria.bp4 is not None
        assert criteria.bp4.applied is True

    def test_protein_domain_gets_pm1(self) -> None:
        """Variant in a protein domain → PM1."""
        annotation = VariantAnnotation(
            ensembl_vep=EnsemblVEPAnnotation(
                found=True, protein_domain="P53_DNA-binding",
            ),
        )
        criteria = _evaluate_criteria_from_evidence(self._make_variant(), annotation)
        assert criteria.pm1 is not None
        assert criteria.pm1.applied is True

    def test_full_pathogenic_evidence(self) -> None:
        """Rare + ClinVar pathogenic + computational deleterious + domain → classification."""
        annotation = VariantAnnotation(
            clinvar=ClinVarAnnotation(
                found=True, clinical_significance="Pathogenic",
                review_stars=3, submitter_count=15,
            ),
            gnomad=GnomADFrequency(found=True, overall_af=0.000005),
            ensembl_vep=EnsemblVEPAnnotation(
                found=True, consequence_type="missense_variant", impact="MODERATE",
                sift_prediction="deleterious", sift_score=0.0,
                polyphen_prediction="probably_damaging", polyphen_score=0.999,
                protein_domain="P53_DNA-binding",
            ),
        )
        criteria = _evaluate_criteria_from_evidence(self._make_variant(), annotation)

        # Should have: PM2, PP5, PP3, PM1 = 2 Moderate + 2 Supporting
        applied = criteria.get_applied_codes()
        codes = [c.code for c in applied]
        assert "PM2" in codes
        assert "PP5" in codes
        assert "PP3" in codes
        assert "PM1" in codes

        result, rule = classify(criteria)
        # 2 Moderate + 2 Supporting = Likely Pathogenic
        assert result.value == "Likely Pathogenic"

    def test_default_annotation_gets_pm2_for_absent_gnomad(self) -> None:
        """Default annotation (gnomAD not found) → PM2 applied.

        This is correct behavior: if gnomAD doesn't find the variant,
        it's absent from population databases.
        """
        annotation = VariantAnnotation()
        criteria = _evaluate_criteria_from_evidence(self._make_variant(), annotation)
        applied = criteria.get_applied_codes()
        assert len(applied) == 1
        assert applied[0].code == "PM2"


class TestCalculateConfidence:
    def test_no_annotation(self) -> None:
        from variantagent.models.classification import ACMGCriteria
        assert _calculate_confidence(None, ACMGCriteria()) == 0.2

    def test_full_data_high_confidence(self) -> None:
        from variantagent.models.classification import ACMGCriteria
        annotation = VariantAnnotation(
            clinvar=ClinVarAnnotation(found=True, review_stars=3),
            gnomad=GnomADFrequency(found=True),
            ensembl_vep=EnsemblVEPAnnotation(found=True),
        )
        confidence = _calculate_confidence(annotation, ACMGCriteria())
        assert confidence >= 0.7
