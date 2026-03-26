"""Variant data models."""

from enum import Enum

from pydantic import BaseModel, Field


class VariantType(str, Enum):
    """Type of genetic variant."""

    SNV = "snv"
    INSERTION = "insertion"
    DELETION = "deletion"
    INDEL = "indel"
    MNV = "mnv"


class Variant(BaseModel):
    """A genetic variant with genomic coordinates."""

    chromosome: str = Field(..., description="Chromosome (e.g., 'chr17', '17')")
    position: int = Field(..., gt=0, description="1-based genomic position")
    reference: str = Field(..., min_length=1, description="Reference allele")
    alternate: str = Field(..., min_length=1, description="Alternate allele")
    gene: str | None = Field(default=None, description="Gene symbol (e.g., 'TP53')")
    variant_type: VariantType | None = Field(default=None, description="Variant type")
    rsid: str | None = Field(default=None, description="dbSNP rsID (e.g., 'rs121913343')")
    hgvs_c: str | None = Field(default=None, description="HGVS coding notation")
    hgvs_p: str | None = Field(default=None, description="HGVS protein notation")
    quality: float | None = Field(default=None, ge=0, description="Variant call quality (QUAL)")
    depth: int | None = Field(default=None, ge=0, description="Read depth at position")
    allele_frequency: float | None = Field(
        default=None, ge=0, le=1, description="Variant allele frequency in the sample"
    )

    @property
    def normalized_chromosome(self) -> str:
        """Return chromosome without 'chr' prefix."""
        return self.chromosome.removeprefix("chr")

    @property
    def variant_id(self) -> str:
        """Return a unique identifier for the variant."""
        return f"{self.chromosome}:{self.position}{self.reference}>{self.alternate}"

    def classify_type(self) -> VariantType:
        """Determine variant type from ref/alt alleles."""
        if len(self.reference) == 1 and len(self.alternate) == 1:
            return VariantType.SNV
        if len(self.reference) < len(self.alternate):
            return VariantType.INSERTION
        if len(self.reference) > len(self.alternate):
            return VariantType.DELETION
        if len(self.reference) > 1 and len(self.reference) == len(self.alternate):
            return VariantType.MNV
        return VariantType.INDEL


class VariantInput(BaseModel):
    """User input for variant analysis — supports VCF file path or manual entry."""

    vcf_path: str | None = Field(default=None, description="Path to VCF file")
    variants: list[Variant] | None = Field(default=None, description="Manually specified variants")
    sample_id: str | None = Field(default=None, description="Sample identifier")
    batch_id: str | None = Field(default=None, description="Batch identifier for batch comparison")

    def get_variants(self) -> list[Variant]:
        """Return variants from either VCF or manual input."""
        if self.variants:
            return self.variants
        if self.vcf_path:
            # VCF parsing is handled by the vcf_parser tool
            raise ValueError("VCF parsing must be done via the vcf_parser tool")
        raise ValueError("Either vcf_path or variants must be provided")
