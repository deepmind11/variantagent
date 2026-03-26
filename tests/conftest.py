"""Shared test fixtures for VariantAgent."""

import pytest

from variantagent.models.variant import Variant, VariantType


@pytest.fixture
def tp53_missense() -> Variant:
    """TP53 R175H — one of the most common cancer mutations."""
    return Variant(
        chromosome="chr17",
        position=7674220,
        reference="G",
        alternate="A",
        gene="TP53",
        variant_type=VariantType.SNV,
        rsid="rs28934578",
        hgvs_p="p.R175H",
        quality=500.0,
        depth=150,
        allele_frequency=0.45,
    )


@pytest.fixture
def brca1_frameshift() -> Variant:
    """BRCA1 frameshift — known pathogenic."""
    return Variant(
        chromosome="chr17",
        position=43091434,
        reference="TG",
        alternate="T",
        gene="BRCA1",
        variant_type=VariantType.DELETION,
        quality=300.0,
        depth=100,
        allele_frequency=0.50,
    )


@pytest.fixture
def benign_snp() -> Variant:
    """A common benign SNP with high population frequency."""
    return Variant(
        chromosome="chr1",
        position=100000,
        reference="A",
        alternate="G",
        gene="UNKNOWN",
        variant_type=VariantType.SNV,
        quality=999.0,
        depth=200,
        allele_frequency=0.50,
    )


@pytest.fixture
def sample_flagstat_text() -> str:
    """Realistic samtools flagstat output for a passing sample."""
    return """10000000 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
50000 + 0 supplementary
500000 + 0 duplicates
500000 + 0 primary duplicates
9800000 + 0 mapped (98.00% : N/A)
9800000 + 0 primary mapped (98.00% : N/A)
9950000 + 0 paired in sequencing
4975000 + 0 read1
4975000 + 0 read2
9500000 + 0 properly paired (95.48% : N/A)
9700000 + 0 with itself and mate mapped
100000 + 0 singletons (1.01% : N/A)
150000 + 0 with mate mapped to a different chr
50000 + 0 with mate mapped to a different chr (mapQ>=5)"""


@pytest.fixture
def failing_flagstat_text() -> str:
    """Samtools flagstat output for a failing sample (high duplication, low mapping)."""
    return """5000000 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
10000 + 0 supplementary
3000000 + 0 duplicates
3000000 + 0 primary duplicates
4000000 + 0 mapped (80.00% : N/A)
4000000 + 0 primary mapped (80.00% : N/A)
4990000 + 0 paired in sequencing
2495000 + 0 read1
2495000 + 0 read2
3500000 + 0 properly paired (70.14% : N/A)
3800000 + 0 with itself and mate mapped
500000 + 0 singletons (10.02% : N/A)
300000 + 0 with mate mapped to a different chr
100000 + 0 with mate mapped to a different chr (mapQ>=5)"""
