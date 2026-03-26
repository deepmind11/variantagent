"""QC Agent: Evaluates sequencing quality around the variant site.

This agent assesses whether the variant call is reliable based on quality metrics.
It has domain-specific thresholds for coverage, duplication, and mapping quality
derived from real production experience with clinical genomics pipelines.

System Prompt Role: Sequencing QC specialist.
Tools: parse_multiqc_json, parse_flagstat, coverage_checker
Distinct Because: Only agent that touches raw QC data. Has domain-specific thresholds.
"""

from __future__ import annotations

from variantagent.models.qc_metrics import QCAssessment, QCIssue, QCStatus

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

# Failure taxonomy: 10+ distinct, realistic QC failure modes
# Each mode has biologically accurate descriptions from production experience
FAILURE_TAXONOMY = {
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


# TODO: Implement QC assessment logic using tools and thresholds above
# The agent should:
# 1. Parse available QC data (MultiQC JSON, flagstat, Picard metrics)
# 2. Compare metrics against thresholds
# 3. Identify the specific failure mode from the taxonomy
# 4. Provide domain-expert level reasoning (not just "coverage is low")
# 5. Return a QCAssessment with issues, causes, and recommendations
