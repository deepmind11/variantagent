"""VCF file parser tool using cyvcf2."""

from __future__ import annotations

from pathlib import Path

from variantagent.models.variant import Variant, VariantType


def parse_vcf(vcf_path: str | Path) -> list[Variant]:
    """Parse a VCF file and return a list of Variant models.

    Args:
        vcf_path: Path to the VCF file (.vcf or .vcf.gz)

    Returns:
        List of parsed Variant objects.

    Raises:
        FileNotFoundError: If VCF file does not exist.
        ValueError: If VCF file is malformed.
    """
    path = Path(vcf_path)
    if not path.exists():
        raise FileNotFoundError(f"VCF file not found: {path}")

    try:
        from cyvcf2 import VCF
    except ImportError as e:
        raise ImportError(
            "cyvcf2 is required for VCF parsing. Install with: pip install cyvcf2"
        ) from e

    variants: list[Variant] = []
    vcf_reader = VCF(str(path))

    for record in vcf_reader:
        for alt in record.ALT:
            variant = Variant(
                chromosome=record.CHROM,
                position=record.POS,
                reference=record.REF,
                alternate=alt,
                quality=record.QUAL,
                depth=record.INFO.get("DP"),
                rsid=record.ID if record.ID and record.ID != "." else None,
            )
            variant.variant_type = variant.classify_type()

            # Extract allele frequency if available
            af_value = record.INFO.get("AF")
            if af_value is not None:
                variant.allele_frequency = float(af_value) if not isinstance(af_value, tuple) else float(af_value[0])

            variants.append(variant)

    vcf_reader.close()
    return variants
