"""Tests for Ensembl VEP client — parsing logic only (no live API calls)."""

from variantagent.models.variant import Variant
from variantagent.tools.ensembl_client import _build_vep_url, _parse_vep_response


class TestBuildVepUrl:
    def test_url_format(self) -> None:
        v = Variant(chromosome="chr17", position=7674220, reference="G", alternate="A")
        url = _build_vep_url(v)
        assert "17:7674220:7674220/A" in url
        assert url.startswith("https://rest.ensembl.org/vep/human/region/")

    def test_strips_chr_prefix(self) -> None:
        v = Variant(chromosome="chr1", position=100, reference="A", alternate="G")
        url = _build_vep_url(v)
        assert "/1:" in url
        assert "chr" not in url.split("region/")[1]


class TestParseVepResponse:
    def test_parse_missense_variant(self) -> None:
        data = [
            {
                "most_severe_consequence": "missense_variant",
                "transcript_consequences": [
                    {
                        "consequence_terms": ["missense_variant"],
                        "impact": "MODERATE",
                        "gene_symbol": "TP53",
                        "gene_id": "ENSG00000141510",
                        "transcript_id": "ENST00000269305",
                        "biotype": "protein_coding",
                        "amino_acids": "R/H",
                        "codons": "cGc/cAc",
                        "sift_prediction": "deleterious",
                        "sift_score": 0.0,
                        "polyphen_prediction": "probably_damaging",
                        "polyphen_score": 0.999,
                        "domains": [
                            {"db": "Pfam", "name": "P53_DNA-binding"},
                        ],
                        "exon": "5/11",
                    }
                ],
            }
        ]
        result = _parse_vep_response(data)
        assert result.found is True
        assert result.consequence_type == "missense_variant"
        assert result.impact == "MODERATE"
        assert result.gene_symbol == "TP53"
        assert result.amino_acid_change == "R/H"
        assert result.sift_prediction == "deleterious"
        assert result.sift_score == 0.0
        assert result.polyphen_prediction == "probably_damaging"
        assert result.polyphen_score == 0.999
        assert result.protein_domain == "P53_DNA-binding"
        assert result.exon == "5/11"

    def test_parse_empty_response(self) -> None:
        result = _parse_vep_response([])
        assert result.found is False

    def test_parse_intergenic_variant(self) -> None:
        data = [
            {
                "most_severe_consequence": "intergenic_variant",
                "intergenic_consequences": [
                    {"consequence_terms": ["intergenic_variant"], "impact": "MODIFIER"}
                ],
            }
        ]
        result = _parse_vep_response(data)
        assert result.found is True
        assert result.consequence_type == "intergenic_variant"
        assert result.impact == "MODIFIER"

    def test_parse_synonymous_no_sift(self) -> None:
        data = [
            {
                "most_severe_consequence": "synonymous_variant",
                "transcript_consequences": [
                    {
                        "consequence_terms": ["synonymous_variant"],
                        "impact": "LOW",
                        "gene_symbol": "BRCA1",
                        "transcript_id": "ENST00000357654",
                        "biotype": "protein_coding",
                    }
                ],
            }
        ]
        result = _parse_vep_response(data)
        assert result.found is True
        assert result.sift_prediction is None
        assert result.polyphen_prediction is None
        assert result.protein_domain is None
