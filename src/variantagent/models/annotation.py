"""Variant annotation data models from public databases."""

from pydantic import BaseModel, Field


class ClinVarAnnotation(BaseModel):
    """ClinVar variant annotation."""

    variation_id: str | None = Field(default=None, description="ClinVar Variation ID")
    clinical_significance: str | None = Field(
        default=None,
        description="Clinical significance (e.g., 'Pathogenic', 'Likely benign', 'VUS')",
    )
    review_status: str | None = Field(
        default=None,
        description="Review status (e.g., 'criteria provided, multiple submitters')",
    )
    review_stars: int | None = Field(
        default=None, ge=0, le=4, description="Review status star rating (0-4)"
    )
    conditions: list[str] = Field(
        default_factory=list, description="Associated conditions/diseases"
    )
    submitter_count: int | None = Field(
        default=None, ge=0, description="Number of submitters"
    )
    last_evaluated: str | None = Field(
        default=None, description="Date of last evaluation"
    )
    found: bool = Field(default=False, description="Whether variant was found in ClinVar")


class GnomADFrequency(BaseModel):
    """gnomAD population allele frequencies."""

    overall_af: float | None = Field(
        default=None, ge=0, le=1, description="Overall allele frequency"
    )
    afr_af: float | None = Field(default=None, ge=0, le=1, description="African/African American")
    amr_af: float | None = Field(default=None, ge=0, le=1, description="Latino/Admixed American")
    asj_af: float | None = Field(default=None, ge=0, le=1, description="Ashkenazi Jewish")
    eas_af: float | None = Field(default=None, ge=0, le=1, description="East Asian")
    fin_af: float | None = Field(default=None, ge=0, le=1, description="Finnish")
    nfe_af: float | None = Field(default=None, ge=0, le=1, description="Non-Finnish European")
    sas_af: float | None = Field(default=None, ge=0, le=1, description="South Asian")
    homozygote_count: int | None = Field(default=None, ge=0, description="Number of homozygotes")
    allele_count: int | None = Field(default=None, ge=0, description="Number of alleles observed")
    allele_number: int | None = Field(default=None, ge=0, description="Total alleles genotyped")
    filtering_status: str | None = Field(
        default=None, description="gnomAD filtering status (PASS, etc.)"
    )
    found: bool = Field(default=False, description="Whether variant was found in gnomAD")


class EnsemblVEPAnnotation(BaseModel):
    """Ensembl Variant Effect Predictor annotation."""

    consequence_type: str | None = Field(
        default=None, description="Most severe consequence (e.g., 'missense_variant')"
    )
    impact: str | None = Field(
        default=None, description="Impact category (HIGH, MODERATE, LOW, MODIFIER)"
    )
    gene_symbol: str | None = Field(default=None, description="Gene symbol")
    gene_id: str | None = Field(default=None, description="Ensembl gene ID")
    transcript_id: str | None = Field(default=None, description="Affected transcript")
    biotype: str | None = Field(default=None, description="Transcript biotype")
    amino_acid_change: str | None = Field(
        default=None, description="Amino acid change (e.g., 'R175H')"
    )
    codon_change: str | None = Field(default=None, description="Codon change")
    sift_prediction: str | None = Field(
        default=None, description="SIFT prediction (tolerated/deleterious)"
    )
    sift_score: float | None = Field(default=None, ge=0, le=1)
    polyphen_prediction: str | None = Field(
        default=None, description="PolyPhen prediction (benign/possibly_damaging/probably_damaging)"
    )
    polyphen_score: float | None = Field(default=None, ge=0, le=1)
    protein_domain: str | None = Field(
        default=None, description="Protein domain affected (e.g., 'p53 DNA-binding domain')"
    )
    exon: str | None = Field(default=None, description="Exon number (e.g., '5/11')")
    found: bool = Field(default=False, description="Whether annotation was retrieved")


class VariantAnnotation(BaseModel):
    """Complete annotation for a variant from all database sources."""

    clinvar: ClinVarAnnotation = Field(default_factory=ClinVarAnnotation)
    gnomad: GnomADFrequency = Field(default_factory=GnomADFrequency)
    ensembl_vep: EnsemblVEPAnnotation = Field(default_factory=EnsemblVEPAnnotation)
    pubmed_references: list[str] = Field(
        default_factory=list, description="Relevant PubMed IDs (PMIDs)"
    )
    annotation_errors: list[str] = Field(
        default_factory=list, description="Errors encountered during annotation"
    )
