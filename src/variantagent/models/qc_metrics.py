"""QC metrics data models for sequencing quality assessment."""

from enum import Enum

from pydantic import BaseModel, Field


class QCStatus(str, Enum):
    """QC assessment status."""

    PASS = "pass"
    WARN = "warn"
    FAIL = "fail"


class FlagstatMetrics(BaseModel):
    """Parsed samtools flagstat output."""

    total_reads: int = Field(..., ge=0)
    mapped_reads: int = Field(..., ge=0)
    mapping_rate: float = Field(..., ge=0, le=1)
    duplicates: int = Field(..., ge=0)
    duplication_rate: float = Field(..., ge=0, le=1)
    paired_reads: int = Field(..., ge=0)
    properly_paired: int = Field(..., ge=0)
    properly_paired_rate: float = Field(..., ge=0, le=1)
    singletons: int = Field(..., ge=0)
    singleton_rate: float = Field(..., ge=0, le=1)


class MultiQCMetrics(BaseModel):
    """Parsed MultiQC general stats."""

    sample_id: str
    total_sequences: int | None = Field(default=None, ge=0)
    percent_gc: float | None = Field(default=None, ge=0, le=100)
    avg_sequence_length: float | None = Field(default=None, ge=0)
    percent_duplicates: float | None = Field(default=None, ge=0, le=100)
    percent_fails: float | None = Field(default=None, ge=0, le=100)
    mean_coverage: float | None = Field(default=None, ge=0)
    median_coverage: float | None = Field(default=None, ge=0)
    percent_bases_above_20x: float | None = Field(default=None, ge=0, le=100)
    percent_bases_above_100x: float | None = Field(default=None, ge=0, le=100)
    insert_size_median: float | None = Field(default=None, ge=0)
    percent_adapter: float | None = Field(default=None, ge=0, le=100)


class QCIssue(BaseModel):
    """A specific QC issue detected."""

    metric: str = Field(..., description="The metric that flagged (e.g., 'duplication_rate')")
    observed_value: float = Field(..., description="The observed value")
    threshold: float = Field(..., description="The threshold that was exceeded")
    severity: QCStatus = Field(..., description="WARN or FAIL")
    description: str = Field(..., description="Human-readable description of the issue")
    likely_causes: list[str] = Field(
        ..., description="Likely root causes based on domain expertise, ranked by probability"
    )
    recommended_action: str = Field(..., description="Recommended action to resolve")


class QCAssessment(BaseModel):
    """Complete QC assessment for a sample or variant region."""

    sample_id: str
    overall_status: QCStatus
    flagstat: FlagstatMetrics | None = None
    multiqc: MultiQCMetrics | None = None
    issues: list[QCIssue] = Field(default_factory=list)
    variant_region_coverage: float | None = Field(
        default=None, ge=0, description="Coverage at the specific variant position"
    )
    variant_region_mapping_quality: float | None = Field(
        default=None, ge=0, description="Average mapping quality at variant position"
    )
    reliable_for_interpretation: bool = Field(
        default=True,
        description="Whether QC supports reliable variant interpretation",
    )
    reasoning: str = Field(
        default="", description="QC Agent's reasoning about the assessment"
    )
