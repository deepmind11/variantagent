"""ACMG/AMP variant classification data models."""

from enum import Enum

from pydantic import BaseModel, Field


class EvidenceStrength(str, Enum):
    """ACMG evidence strength levels."""

    VERY_STRONG = "very_strong"
    STRONG = "strong"
    MODERATE = "moderate"
    SUPPORTING = "supporting"


class EvidenceDirection(str, Enum):
    """Whether evidence supports pathogenicity or benign classification."""

    PATHOGENIC = "pathogenic"
    BENIGN = "benign"


class EvidenceCode(BaseModel):
    """A single ACMG evidence criterion assessment."""

    code: str = Field(
        ..., description="ACMG evidence code (e.g., 'PVS1', 'PS1', 'PM2', 'BA1')"
    )
    name: str = Field(..., description="Human-readable name of the criterion")
    direction: EvidenceDirection = Field(..., description="Pathogenic or benign evidence")
    strength: EvidenceStrength = Field(..., description="Evidence strength level")
    applied: bool = Field(..., description="Whether this criterion is met")
    reasoning: str = Field(..., description="Chain-of-thought reasoning for this criterion")
    data_source: str = Field(..., description="Which database/tool provided the evidence")
    confidence: float = Field(
        ..., ge=0, le=1, description="Confidence in this criterion assessment"
    )


class ACMGClassificationResult(str, Enum):
    """ACMG/AMP five-tier classification."""

    PATHOGENIC = "Pathogenic"
    LIKELY_PATHOGENIC = "Likely Pathogenic"
    VUS = "Uncertain Significance"
    LIKELY_BENIGN = "Likely Benign"
    BENIGN = "Benign"


class ACMGCriteria(BaseModel):
    """Complete set of ACMG criteria evaluated for a variant."""

    # Pathogenic criteria
    pvs1: EvidenceCode | None = Field(default=None, description="Null variant in gene where LOF is known mechanism")
    ps1: EvidenceCode | None = Field(default=None, description="Same amino acid change as established pathogenic")
    ps3: EvidenceCode | None = Field(default=None, description="Functional studies supportive")
    pm1: EvidenceCode | None = Field(default=None, description="In mutational hot spot / functional domain")
    pm2: EvidenceCode | None = Field(default=None, description="Absent from population databases")
    pm4: EvidenceCode | None = Field(default=None, description="Protein length change from in-frame indel")
    pm5: EvidenceCode | None = Field(default=None, description="Novel missense at position with known pathogenic")
    pp2: EvidenceCode | None = Field(default=None, description="Missense in gene with low rate of benign missense")
    pp3: EvidenceCode | None = Field(default=None, description="Computational evidence supports deleterious")
    pp5: EvidenceCode | None = Field(default=None, description="Reputable source reports pathogenic")

    # Benign criteria
    ba1: EvidenceCode | None = Field(default=None, description="Allele frequency > 5% (standalone benign)")
    bs1: EvidenceCode | None = Field(default=None, description="Allele frequency greater than expected for disorder")
    bs2: EvidenceCode | None = Field(default=None, description="Observed in healthy adult (dominant) or homozygous (recessive)")
    bp1: EvidenceCode | None = Field(default=None, description="Missense in gene where truncating is mechanism")
    bp4: EvidenceCode | None = Field(default=None, description="Computational evidence supports benign")
    bp6: EvidenceCode | None = Field(default=None, description="Reputable source reports benign")
    bp7: EvidenceCode | None = Field(default=None, description="Silent variant with no splicing impact")

    def get_applied_codes(self) -> list[EvidenceCode]:
        """Return all criteria that were applied (met)."""
        codes = []
        for field_name in self.model_fields:
            value = getattr(self, field_name)
            if isinstance(value, EvidenceCode) and value.applied:
                codes.append(value)
        return codes

    def get_pathogenic_codes(self) -> list[EvidenceCode]:
        """Return applied pathogenic evidence codes."""
        return [c for c in self.get_applied_codes() if c.direction == EvidenceDirection.PATHOGENIC]

    def get_benign_codes(self) -> list[EvidenceCode]:
        """Return applied benign evidence codes."""
        return [c for c in self.get_applied_codes() if c.direction == EvidenceDirection.BENIGN]


class ACMGClassification(BaseModel):
    """Final ACMG/AMP classification for a variant."""

    classification: ACMGClassificationResult
    criteria: ACMGCriteria
    confidence: float = Field(..., ge=0, le=1, description="Overall confidence in classification")
    reasoning: str = Field(..., description="Summary reasoning for the classification")
    applied_codes_summary: list[str] = Field(
        default_factory=list, description="List of applied evidence codes (e.g., ['PM2', 'PP3'])"
    )
    classification_rule: str = Field(
        ...,
        description="The ACMG combining rule that produced this classification "
        "(e.g., '1 Strong + 1 Moderate = Likely Pathogenic')",
    )
