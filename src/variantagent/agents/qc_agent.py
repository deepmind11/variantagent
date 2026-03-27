"""QC Agent: Evaluates sequencing quality around the variant site.

This agent assesses whether the variant call is reliable based on quality metrics.
It has domain-specific thresholds for coverage, duplication, and mapping quality
derived from real production experience with clinical genomics pipelines.

System Prompt Role: Sequencing QC specialist.
Tools: parse_multiqc_json, parse_flagstat, coverage_checker
Distinct Because: Only agent that touches raw QC data. Has domain-specific thresholds.
"""

from __future__ import annotations

import logging

from variantagent.models.qc_metrics import (
    FlagstatMetrics,
    MultiQCMetrics,
    QCAssessment,
    QCIssue,
    QCStatus,
)

logger = logging.getLogger(__name__)

# Domain-specific thresholds from clinical genomics production experience
# These encode knowledge from triaging thousands of real samples
QC_THRESHOLDS = {
    "min_coverage": {
        "fail": 20,
        "warn": 50,
        "description": "Minimum mean coverage depth",
    },
    "max_duplication_rate": {
        "fail": 0.50,
        "warn": 0.30,
        "description": "Maximum duplication rate",
    },
    "min_mapping_rate": {
        "fail": 0.90,
        "warn": 0.95,
        "description": "Minimum mapping rate",
    },
    "min_properly_paired_rate": {
        "fail": 0.85,
        "warn": 0.90,
        "description": "Minimum properly paired rate",
    },
    "max_singleton_rate": {
        "fail": 0.10,
        "warn": 0.05,
        "description": "Maximum singleton rate",
    },
    "min_variant_position_coverage": {
        "fail": 10,
        "warn": 30,
        "description": "Minimum coverage at the specific variant position",
    },
}

# Failure taxonomy: 11 distinct, realistic QC failure modes
# Each mode has biologically accurate descriptions from production experience
FAILURE_TAXONOMY: dict[str, dict[str, str | list[str]]] = {
    "low_coverage_global": {
        "description": "Global low coverage across the panel/genome",
        "likely_causes": [
            "Low DNA input quantity",
            "Underclustering during sequencing",
            "Failed library amplification",
        ],
        "recommended_action": "Check input DNA concentration. If adequate, re-sequence.",
    },
    "low_coverage_regional": {
        "description": "Low coverage at the specific variant position while global coverage is adequate",
        "likely_causes": [
            "GC bias in capture — region has extreme GC content",
            "Capture probe failure for this target region",
            "Repetitive sequence causing mapping ambiguity",
        ],
        "recommended_action": "Check GC content at position. Review capture panel design for this region.",
    },
    "high_duplication": {
        "description": "Elevated duplicate rate indicating low library complexity",
        "likely_causes": [
            "Low DNA input requiring excessive PCR amplification",
            "Library prep protocol issue (over-amplification)",
            "Sample degradation (FFPE, cfDNA with very low input)",
        ],
        "recommended_action": "Check input DNA amount and library complexity metrics. Consider re-extraction.",
    },
    "high_duplication_batch": {
        "description": "High duplication affecting multiple samples in the same batch",
        "likely_causes": [
            "Library prep batch effect — reagent issue or protocol deviation",
            "Systematic over-amplification across the batch",
        ],
        "recommended_action": "Compare with other samples in batch. If batch-wide, flag library prep for review.",
    },
    "contamination": {
        "description": "Evidence of sample contamination from a second individual",
        "likely_causes": [
            "Sample swap during accessioning",
            "Cross-contamination during library prep",
            "Index hopping in multiplexed sequencing",
        ],
        "recommended_action": "Run contamination estimation (VerifyBamID). If confirmed, re-extract from backup.",
    },
    "low_mapping_rate": {
        "description": "Abnormally low percentage of reads mapping to reference",
        "likely_causes": [
            "Contamination with non-human DNA (microbial, adapter dimers)",
            "Wrong reference genome used for alignment",
            "Sample swap with a different species",
        ],
        "recommended_action": "Run BLAST on unmapped reads. Check species of origin.",
    },
    "insert_size_anomaly": {
        "description": "Insert size distribution outside expected range",
        "likely_causes": [
            "DNA degradation (FFPE samples have characteristically short fragments)",
            "Incomplete shearing during library prep",
            "cfDNA from liquid biopsy (expected ~167bp mononucleosomal peak)",
        ],
        "recommended_action": "Review insert size histogram. For FFPE/cfDNA, thresholds differ from standard.",
    },
    "adapter_contamination": {
        "description": "High percentage of reads containing adapter sequences",
        "likely_causes": [
            "Short library fragments where read-through occurs",
            "Adapter dimer contamination from library prep",
        ],
        "recommended_action": "Verify adapter trimming was applied. Check library fragment size distribution.",
    },
    "strand_bias": {
        "description": "Variant supported predominantly by reads from one strand",
        "likely_causes": [
            "Sequencing artifact (common with certain error modes)",
            "DNA damage artifact (oxidative damage creates G>T on one strand)",
            "True biological variant in a difficult region",
        ],
        "recommended_action": "Check strand bias metrics in VCF. If SB > threshold, flag variant call as unreliable.",
    },
    "sample_swap": {
        "description": "Genotype at known SNP positions does not match expected sample identity",
        "likely_causes": [
            "Sample mislabeling during accessioning",
            "Plate position swap during library prep",
            "Barcode collision or index assignment error",
        ],
        "recommended_action": "Run fingerprinting check against expected genotypes. Cross-reference with sample manifest.",
    },
    "allelic_imbalance": {
        "description": "Unexpected allele frequency for the variant call",
        "likely_causes": [
            "Subclonal variant (somatic, expected in oncology)",
            "Copy number alteration at the locus",
            "Mosaic variant",
            "Contamination at low level",
        ],
        "recommended_action": "For oncology: check if VAF is consistent with tumor purity. For germline: investigate CNV.",
    },
}


def assess_flagstat(flagstat: FlagstatMetrics) -> list[QCIssue]:
    """Evaluate flagstat metrics against thresholds and return issues found.

    This is the core QC logic — it encodes domain expertise about what
    metric values indicate which failure modes.
    """
    issues: list[QCIssue] = []

    # Check mapping rate
    if flagstat.mapping_rate < QC_THRESHOLDS["min_mapping_rate"]["fail"]:
        taxonomy = FAILURE_TAXONOMY["low_mapping_rate"]
        issues.append(
            QCIssue(
                metric="mapping_rate",
                observed_value=flagstat.mapping_rate,
                threshold=QC_THRESHOLDS["min_mapping_rate"]["fail"],
                severity=QCStatus.FAIL,
                description=str(taxonomy["description"]),
                likely_cause=str(taxonomy["likely_causes"][0]) if isinstance(taxonomy["likely_causes"], list) else str(taxonomy["likely_causes"]),
                recommended_action=str(taxonomy["recommended_action"]),
            )
        )
    elif flagstat.mapping_rate < QC_THRESHOLDS["min_mapping_rate"]["warn"]:
        issues.append(
            QCIssue(
                metric="mapping_rate",
                observed_value=flagstat.mapping_rate,
                threshold=QC_THRESHOLDS["min_mapping_rate"]["warn"],
                severity=QCStatus.WARN,
                description="Mapping rate below optimal threshold",
                likely_cause="Possible low-level contamination or difficult library",
                recommended_action="Monitor — acceptable but investigate if other metrics also flag.",
            )
        )

    # Check duplication rate
    if flagstat.duplication_rate > QC_THRESHOLDS["max_duplication_rate"]["fail"]:
        taxonomy = FAILURE_TAXONOMY["high_duplication"]
        issues.append(
            QCIssue(
                metric="duplication_rate",
                observed_value=flagstat.duplication_rate,
                threshold=QC_THRESHOLDS["max_duplication_rate"]["fail"],
                severity=QCStatus.FAIL,
                description=str(taxonomy["description"]),
                likely_cause=str(taxonomy["likely_causes"][0]) if isinstance(taxonomy["likely_causes"], list) else str(taxonomy["likely_causes"]),
                recommended_action=str(taxonomy["recommended_action"]),
            )
        )
    elif flagstat.duplication_rate > QC_THRESHOLDS["max_duplication_rate"]["warn"]:
        issues.append(
            QCIssue(
                metric="duplication_rate",
                observed_value=flagstat.duplication_rate,
                threshold=QC_THRESHOLDS["max_duplication_rate"]["warn"],
                severity=QCStatus.WARN,
                description="Elevated duplicate rate — library complexity may be low",
                likely_cause="Moderate PCR over-amplification or borderline input DNA",
                recommended_action="Note for interpretation — effective coverage is reduced.",
            )
        )

    # Check properly paired rate
    if flagstat.properly_paired_rate < QC_THRESHOLDS["min_properly_paired_rate"]["fail"]:
        issues.append(
            QCIssue(
                metric="properly_paired_rate",
                observed_value=flagstat.properly_paired_rate,
                threshold=QC_THRESHOLDS["min_properly_paired_rate"]["fail"],
                severity=QCStatus.FAIL,
                description="Low properly-paired rate indicates structural issues",
                likely_cause="Possible chimeric reads, contamination, or alignment artifacts",
                recommended_action="Investigate insert size distribution and check for chimeras.",
            )
        )
    elif flagstat.properly_paired_rate < QC_THRESHOLDS["min_properly_paired_rate"]["warn"]:
        issues.append(
            QCIssue(
                metric="properly_paired_rate",
                observed_value=flagstat.properly_paired_rate,
                threshold=QC_THRESHOLDS["min_properly_paired_rate"]["warn"],
                severity=QCStatus.WARN,
                description="Properly-paired rate below optimal threshold",
                likely_cause="Minor alignment issues or library prep anomaly",
                recommended_action="Monitor — acceptable but note for interpretation.",
            )
        )

    # Check singleton rate
    if flagstat.singleton_rate > QC_THRESHOLDS["max_singleton_rate"]["fail"]:
        issues.append(
            QCIssue(
                metric="singleton_rate",
                observed_value=flagstat.singleton_rate,
                threshold=QC_THRESHOLDS["max_singleton_rate"]["fail"],
                severity=QCStatus.FAIL,
                description="High singleton rate indicates mate-pair issues",
                likely_cause="Library prep failure, chimeric reads, or contamination",
                recommended_action="Check insert size and mate mapping. Consider re-prep.",
            )
        )
    elif flagstat.singleton_rate > QC_THRESHOLDS["max_singleton_rate"]["warn"]:
        issues.append(
            QCIssue(
                metric="singleton_rate",
                observed_value=flagstat.singleton_rate,
                threshold=QC_THRESHOLDS["max_singleton_rate"]["warn"],
                severity=QCStatus.WARN,
                description="Elevated singleton rate",
                likely_cause="Minor mate-pair mapping issues",
                recommended_action="Monitor — typically not actionable alone.",
            )
        )

    return issues


def assess_multiqc(multiqc: MultiQCMetrics) -> list[QCIssue]:
    """Evaluate MultiQC metrics against thresholds."""
    issues: list[QCIssue] = []

    # Check coverage
    coverage = multiqc.mean_coverage
    if coverage is not None:
        if coverage < QC_THRESHOLDS["min_coverage"]["fail"]:
            taxonomy = FAILURE_TAXONOMY["low_coverage_global"]
            issues.append(
                QCIssue(
                    metric="mean_coverage",
                    observed_value=coverage,
                    threshold=QC_THRESHOLDS["min_coverage"]["fail"],
                    severity=QCStatus.FAIL,
                    description=str(taxonomy["description"]),
                    likely_cause=str(taxonomy["likely_causes"][0]) if isinstance(taxonomy["likely_causes"], list) else str(taxonomy["likely_causes"]),
                    recommended_action=str(taxonomy["recommended_action"]),
                )
            )
        elif coverage < QC_THRESHOLDS["min_coverage"]["warn"]:
            issues.append(
                QCIssue(
                    metric="mean_coverage",
                    observed_value=coverage,
                    threshold=QC_THRESHOLDS["min_coverage"]["warn"],
                    severity=QCStatus.WARN,
                    description="Coverage below optimal threshold",
                    likely_cause="Borderline DNA input or minor sequencing underperformance",
                    recommended_action="Variant calls in low-complexity regions may be less reliable.",
                )
            )

    # Check adapter contamination
    if multiqc.percent_adapter is not None and multiqc.percent_adapter > 5.0:
        taxonomy = FAILURE_TAXONOMY["adapter_contamination"]
        severity = QCStatus.FAIL if multiqc.percent_adapter > 20.0 else QCStatus.WARN
        issues.append(
            QCIssue(
                metric="percent_adapter",
                observed_value=multiqc.percent_adapter,
                threshold=5.0,
                severity=severity,
                description=str(taxonomy["description"]),
                likely_cause=str(taxonomy["likely_causes"][0]) if isinstance(taxonomy["likely_causes"], list) else str(taxonomy["likely_causes"]),
                recommended_action=str(taxonomy["recommended_action"]),
            )
        )

    return issues


def run_qc_assessment(
    sample_id: str,
    flagstat: FlagstatMetrics | None = None,
    multiqc: MultiQCMetrics | None = None,
    variant_region_coverage: float | None = None,
) -> QCAssessment:
    """Run full QC assessment on a sample.

    This is the main entry point for the QC Agent. It evaluates all available
    QC data and produces a structured assessment.

    Args:
        sample_id: Sample identifier.
        flagstat: Parsed flagstat metrics (optional).
        multiqc: Parsed MultiQC metrics (optional).
        variant_region_coverage: Coverage at the specific variant position (optional).

    Returns:
        Complete QC assessment with issues, status, and recommendations.
    """
    all_issues: list[QCIssue] = []

    if flagstat is not None:
        all_issues.extend(assess_flagstat(flagstat))

    if multiqc is not None:
        all_issues.extend(assess_multiqc(multiqc))

    # Check variant-position-specific coverage
    if variant_region_coverage is not None:
        if variant_region_coverage < QC_THRESHOLDS["min_variant_position_coverage"]["fail"]:
            taxonomy = FAILURE_TAXONOMY["low_coverage_regional"]
            all_issues.append(
                QCIssue(
                    metric="variant_position_coverage",
                    observed_value=variant_region_coverage,
                    threshold=QC_THRESHOLDS["min_variant_position_coverage"]["fail"],
                    severity=QCStatus.FAIL,
                    description=str(taxonomy["description"]),
                    likely_cause=str(taxonomy["likely_causes"][0]) if isinstance(taxonomy["likely_causes"], list) else str(taxonomy["likely_causes"]),
                    recommended_action=str(taxonomy["recommended_action"]),
                )
            )
        elif variant_region_coverage < QC_THRESHOLDS["min_variant_position_coverage"]["warn"]:
            all_issues.append(
                QCIssue(
                    metric="variant_position_coverage",
                    observed_value=variant_region_coverage,
                    threshold=QC_THRESHOLDS["min_variant_position_coverage"]["warn"],
                    severity=QCStatus.WARN,
                    description="Coverage at variant position is below optimal",
                    likely_cause="Regional capture or GC bias effect",
                    recommended_action="Interpret variant with caution — low coverage reduces call confidence.",
                )
            )

    # Determine overall status
    has_fail = any(issue.severity == QCStatus.FAIL for issue in all_issues)
    has_warn = any(issue.severity == QCStatus.WARN for issue in all_issues)

    if has_fail:
        overall_status = QCStatus.FAIL
    elif has_warn:
        overall_status = QCStatus.WARN
    else:
        overall_status = QCStatus.PASS

    # Determine if QC supports reliable interpretation
    reliable = not has_fail
    if variant_region_coverage is not None and variant_region_coverage < QC_THRESHOLDS["min_variant_position_coverage"]["fail"]:
        reliable = False

    # Build reasoning summary
    if not all_issues:
        reasoning = "All QC metrics within acceptable thresholds. Variant call is reliable."
    else:
        fail_count = sum(1 for i in all_issues if i.severity == QCStatus.FAIL)
        warn_count = sum(1 for i in all_issues if i.severity == QCStatus.WARN)
        parts = []
        if fail_count:
            parts.append(f"{fail_count} FAIL")
        if warn_count:
            parts.append(f"{warn_count} WARN")
        reasoning = f"QC assessment: {', '.join(parts)} issues detected. "
        if not reliable:
            reasoning += "Variant interpretation may not be reliable due to QC failures."
        else:
            reasoning += "Issues are warnings only — interpretation can proceed with caution."

    logger.info("QC assessment for %s: %s (%d issues)", sample_id, overall_status.value, len(all_issues))

    return QCAssessment(
        sample_id=sample_id,
        overall_status=overall_status,
        flagstat=flagstat,
        multiqc=multiqc,
        issues=all_issues,
        variant_region_coverage=variant_region_coverage,
        reliable_for_interpretation=reliable,
        reasoning=reasoning,
    )
