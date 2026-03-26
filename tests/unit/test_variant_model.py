"""Tests for Variant data models."""

import pytest

from variantagent.models.variant import Variant, VariantType


class TestVariant:
    def test_snv_classification(self) -> None:
        v = Variant(chromosome="chr1", position=100, reference="A", alternate="G")
        assert v.classify_type() == VariantType.SNV

    def test_insertion_classification(self) -> None:
        v = Variant(chromosome="chr1", position=100, reference="A", alternate="ATG")
        assert v.classify_type() == VariantType.INSERTION

    def test_deletion_classification(self) -> None:
        v = Variant(chromosome="chr1", position=100, reference="ATG", alternate="A")
        assert v.classify_type() == VariantType.DELETION

    def test_mnv_classification(self) -> None:
        v = Variant(chromosome="chr1", position=100, reference="AT", alternate="GC")
        assert v.classify_type() == VariantType.MNV

    def test_normalized_chromosome(self) -> None:
        v = Variant(chromosome="chr17", position=100, reference="A", alternate="G")
        assert v.normalized_chromosome == "17"

    def test_normalized_chromosome_no_prefix(self) -> None:
        v = Variant(chromosome="17", position=100, reference="A", alternate="G")
        assert v.normalized_chromosome == "17"

    def test_variant_id(self) -> None:
        v = Variant(chromosome="chr17", position=7674220, reference="G", alternate="A")
        assert v.variant_id == "chr17:7674220G>A"

    def test_validation_rejects_zero_position(self) -> None:
        with pytest.raises(ValueError):
            Variant(chromosome="chr1", position=0, reference="A", alternate="G")

    def test_validation_rejects_empty_allele(self) -> None:
        with pytest.raises(ValueError):
            Variant(chromosome="chr1", position=100, reference="", alternate="G")

    def test_allele_frequency_bounds(self) -> None:
        with pytest.raises(ValueError):
            Variant(
                chromosome="chr1",
                position=100,
                reference="A",
                alternate="G",
                allele_frequency=1.5,
            )
