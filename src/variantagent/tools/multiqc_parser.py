"""Parser for MultiQC general stats JSON output."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from variantagent.models.qc_metrics import MultiQCMetrics


def parse_multiqc_json(multiqc_path: str | Path) -> list[MultiQCMetrics]:
    """Parse MultiQC general stats JSON into structured metrics.

    Handles the multiqc_data/multiqc_general_stats.txt or
    multiqc_data.json format.

    Args:
        multiqc_path: Path to the MultiQC JSON data file.

    Returns:
        List of MultiQCMetrics, one per sample.

    Raises:
        FileNotFoundError: If file does not exist.
        ValueError: If file format is unexpected.
    """
    path = Path(multiqc_path)
    if not path.exists():
        raise FileNotFoundError(f"MultiQC file not found: {path}")

    content = path.read_text()
    data = json.loads(content)

    return parse_multiqc_data(data)


def parse_multiqc_data(data: dict[str, Any]) -> list[MultiQCMetrics]:
    """Parse MultiQC data from a parsed JSON dict.

    Supports both the top-level multiqc_data.json format and the
    report_general_stats format.

    Args:
        data: Parsed MultiQC JSON data.

    Returns:
        List of MultiQCMetrics, one per sample.
    """
    metrics_list: list[MultiQCMetrics] = []

    # Handle multiqc_data.json format (report_general_stats)
    general_stats = data.get("report_general_stats_data", data.get("general_stats", []))

    if isinstance(general_stats, list):
        for stats_block in general_stats:
            if isinstance(stats_block, dict):
                for sample_id, stats in stats_block.items():
                    metrics_list.append(_extract_metrics(sample_id, stats))
    elif isinstance(general_stats, dict):
        for sample_id, stats in general_stats.items():
            metrics_list.append(_extract_metrics(sample_id, stats))

    return metrics_list


def _extract_metrics(sample_id: str, stats: dict[str, Any]) -> MultiQCMetrics:
    """Extract metrics from a single sample's stats dict."""
    return MultiQCMetrics(
        sample_id=sample_id,
        total_sequences=_safe_int(stats.get("total_sequences")),
        percent_gc=_safe_float(stats.get("percent_gc")),
        avg_sequence_length=_safe_float(stats.get("avg_sequence_length")),
        percent_duplicates=_safe_float(stats.get("percent_duplicates")),
        percent_fails=_safe_float(stats.get("percent_fails")),
        mean_coverage=_safe_float(
            stats.get("mean_coverage", stats.get("mosdepth_mean_coverage"))
        ),
        median_coverage=_safe_float(
            stats.get("median_coverage", stats.get("mosdepth_median_coverage"))
        ),
        percent_bases_above_20x=_safe_float(stats.get("20_x_pc")),
        percent_bases_above_100x=_safe_float(stats.get("100_x_pc")),
        insert_size_median=_safe_float(
            stats.get("summed_median", stats.get("median_insert_size"))
        ),
        percent_adapter=_safe_float(stats.get("percent_adapter")),
    )


def _safe_float(value: Any) -> float | None:
    """Safely convert a value to float."""
    if value is None:
        return None
    try:
        return float(value)
    except (ValueError, TypeError):
        return None


def _safe_int(value: Any) -> int | None:
    """Safely convert a value to int."""
    if value is None:
        return None
    try:
        return int(value)
    except (ValueError, TypeError):
        return None
