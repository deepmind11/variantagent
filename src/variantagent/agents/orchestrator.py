"""Orchestrator Agent: LangGraph state machine that plans and routes analysis."""

from __future__ import annotations

import uuid
from typing import TypedDict

from variantagent.models.annotation import VariantAnnotation
from variantagent.models.classification import ACMGClassification
from variantagent.models.qc_metrics import QCAssessment
from variantagent.models.report import TriageReport
from variantagent.models.variant import Variant


class AnalysisState(TypedDict):
    """LangGraph state for the variant analysis workflow."""

    trace_id: str
    variant: Variant
    sample_id: str | None
    batch_id: str | None
    plan: list[str]
    qc_assessment: QCAssessment | None
    annotation: VariantAnnotation | None
    classification: ACMGClassification | None
    report: TriageReport | None
    errors: list[str]
    current_step: int


def create_initial_state(variant: Variant, sample_id: str | None = None, batch_id: str | None = None) -> AnalysisState:
    """Create initial state for a new analysis run."""
    return AnalysisState(
        trace_id=str(uuid.uuid4()),
        variant=variant,
        sample_id=sample_id,
        batch_id=batch_id,
        plan=[],
        qc_assessment=None,
        annotation=None,
        classification=None,
        report=None,
        errors=[],
        current_step=0,
    )


# TODO: Implement LangGraph StateGraph with:
# 1. plan_node: Creates analysis plan based on variant type and available data
# 2. qc_node: Runs QC Agent (skipped if no QC data available)
# 3. route_after_qc: Dynamic routing — if QC fails, skip to report with warning
# 4. annotation_node: Runs Annotation Agent (ClinVar, gnomAD, Ensembl VEP)
# 5. route_after_annotation: If novel variant found, trigger Literature Agent
# 6. literature_node: Runs Literature Agent (PubMed + RAG over ACMG guidelines)
# 7. classification_node: Runs Classification Agent (ACMG rule engine)
# 8. review_node: Runs Reviewer Agent (self-evaluation, contradiction detection)
# 9. hitl_node: Human-in-the-loop checkpoint (confidence-gated)
# 10. report_node: Generates final TriageReport with provenance
