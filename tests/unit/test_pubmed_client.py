"""Tests for PubMed client pure functions."""

from variantagent.tools.pubmed_client import PubMedArticle, _build_search_queries


class TestPubMedArticle:
    def test_citation_with_authors(self) -> None:
        article = PubMedArticle(
            pmid="12345678",
            title="BRCA1 variant study",
            journal="Nature Genetics",
            year="2023",
            authors=["Smith J", "Doe A"],
        )
        citation = article.citation()
        assert "Smith J et al." in citation
        assert "2023" in citation
        assert "BRCA1 variant study" in citation
        assert "Nature Genetics" in citation
        assert "PMID: 12345678" in citation

    def test_citation_no_authors(self) -> None:
        article = PubMedArticle(
            pmid="99999",
            title="Some title",
        )
        citation = article.citation()
        assert "Unknown et al." in citation

    def test_default_attributes(self) -> None:
        article = PubMedArticle(pmid="111", title="Test")
        assert article.journal == ""
        assert article.year == ""
        assert article.authors == []

    def test_authors_default_none_becomes_empty_list(self) -> None:
        article = PubMedArticle(pmid="222", title="Test", authors=None)
        assert article.authors == []


class TestBuildSearchQueries:
    def test_with_gene_and_hgvs_p(self) -> None:
        queries = _build_search_queries("TP53", "p.R175H", "chr17:7674220:G:A")
        assert len(queries) >= 1
        # First query should be most specific (p. prefix stripped)
        assert "TP53" in queries[0]
        assert "R175H" in queries[0]

    def test_with_gene_only(self) -> None:
        queries = _build_search_queries("BRCA1", None, "chr17:43091434:TG:T")
        assert len(queries) >= 1
        for q in queries:
            assert "BRCA1" in q

    def test_returns_list_of_strings(self) -> None:
        queries = _build_search_queries("EGFR", "p.L858R", "chr7:55259515:T:G")
        assert isinstance(queries, list)
        assert all(isinstance(q, str) for q in queries)

    def test_progressive_broadening(self) -> None:
        queries = _build_search_queries("KRAS", "p.G12D", "chr12:25398281:C:T")
        # Should have multiple queries (specific → broad)
        assert len(queries) >= 2
