"""Report and provenance data models."""

from datetime import datetime

from pydantic import BaseModel, Field

from variantagent.models.annotation import VariantAnnotation
from variantagent.models.classification import ACMGClassification
from variantagent.models.qc_metrics import QCAssessment
from variantagent.models.variant import Variant


class ProvenanceEntry(BaseModel):
    """A single step in the provenance trail — traces every conclusion to its source."""

    step: int = Field(..., description="Step number in the analysis")
    agent: str = Field(..., description="Which agent performed this step")
    action: str = Field(..., description="What the agent did")
    tool_called: str | None = Field(default=None, description="Tool/API called")
    input_summary: str = Field(..., description="Summary of input to the tool")
    output_summary: str = Field(..., description="Summary of output from the tool")
    data_source: str | None = Field(default=None, description="External data source queried")
    timestamp: datetime = Field(default_factory=datetime.now)
    duration_ms: int | None = Field(default=None, ge=0, description="Duration of this step in ms")
    tokens_used: int | None = Field(default=None, ge=0, description="LLM tokens consumed")
    error: str | None = Field(default=None, description="Error message if step failed")
    recovered: bool = Field(
        default=False, description="Whether error recovery was attempted and succeeded"
    )


class ReviewerFinding(BaseModel):
    """A finding from the Reviewer Agent's self-evaluation."""

    claim: str = Field(..., description="The claim being evaluated")
    supported: bool = Field(..., description="Whether the claim is supported by evidence")
    source_references: list[str] = Field(
        default_factory=list, description="Sources that support or contradict the claim"
    )
    concern: str | None = Field(
        default=None, description="Concern raised by the reviewer (if any)"
    )
    hallucination_risk: str = Field(
        default="low", description="Hallucination risk: low, medium, high"
    )


class TriageReport(BaseModel):
    """Complete variant interpretation report with full provenance."""

    # Identifiers
    trace_id: str = Field(..., description="Unique trace ID for this analysis run")
    timestamp: datetime = Field(default_factory=datetime.now)

    # Input
    variant: Variant
    sample_id: str | None = None
    batch_id: str | None = None

    # Analysis results
    qc_assessment: QCAssessment | None = None
    annotation: VariantAnnotation | None = None
    classification: ACMGClassification | None = None

    # Self-evaluation
    reviewer_findings: list[ReviewerFinding] = Field(default_factory=list)
    overall_confidence: float = Field(
        ..., ge=0, le=1, description="Overall confidence in the report"
    )
    requires_human_review: bool = Field(
        default=False, description="Whether this report was flagged for human review"
    )
    human_review_reason: str | None = Field(
        default=None, description="Why human review is needed"
    )

    # Provenance
    provenance: list[ProvenanceEntry] = Field(
        default_factory=list, description="Complete audit trail of every step"
    )
    analysis_plan: list[str] = Field(
        default_factory=list,
        description="The plan the Orchestrator created before execution",
    )

    # Cost and performance
    total_duration_ms: int | None = Field(default=None, ge=0)
    total_tokens_used: int | None = Field(default=None, ge=0)
    estimated_cost_usd: float | None = Field(default=None, ge=0)

    # Summary
    natural_language_summary: str = Field(
        default="", description="Human-readable summary of the findings"
    )
    limitations: list[str] = Field(
        default_factory=list,
        description="Honest limitations of this analysis",
    )
