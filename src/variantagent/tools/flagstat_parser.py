"""Parser for samtools flagstat output."""

from __future__ import annotations

import re
from pathlib import Path

from variantagent.models.qc_metrics import FlagstatMetrics


def parse_flagstat(flagstat_path: str | Path) -> FlagstatMetrics:
    """Parse samtools flagstat output into structured metrics.

    Handles the standard samtools flagstat format:
        12345 + 0 in total (QC-passed reads + QC-failed reads)
        100 + 0 secondary
        0 + 0 supplementary
        500 + 0 duplicates
        ...

    Args:
        flagstat_path: Path to the flagstat output file.

    Returns:
        Parsed FlagstatMetrics.

    Raises:
        FileNotFoundError: If file does not exist.
        ValueError: If file format is unexpected.
    """
    path = Path(flagstat_path)
    if not path.exists():
        raise FileNotFoundError(f"Flagstat file not found: {path}")

    content = path.read_text()
    return parse_flagstat_text(content)


def parse_flagstat_text(text: str) -> FlagstatMetrics:
    """Parse flagstat from raw text content.

    Args:
        text: Raw samtools flagstat output text.

    Returns:
        Parsed FlagstatMetrics.
    """
    lines = text.strip().split("\n")

    def extract_count(pattern: str) -> int:
        for line in lines:
            if pattern in line:
                match = re.match(r"(\d+)\s*\+\s*(\d+)", line)
                if match:
                    return int(match.group(1))
        return 0

    total = extract_count("in total")
    duplicates = extract_count("duplicates")
    mapped = extract_count("mapped (")
    paired = extract_count("paired in sequencing")
    properly_paired = extract_count("properly paired")
    singletons = extract_count("singletons")

    return FlagstatMetrics(
        total_reads=total,
        mapped_reads=mapped,
        mapping_rate=mapped / total if total > 0 else 0.0,
        duplicates=duplicates,
        duplication_rate=duplicates / total if total > 0 else 0.0,
        paired_reads=paired,
        properly_paired=properly_paired,
        properly_paired_rate=properly_paired / paired if paired > 0 else 0.0,
        singletons=singletons,
        singleton_rate=singletons / total if total > 0 else 0.0,
    )
