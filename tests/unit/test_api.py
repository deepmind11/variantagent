"""Tests for the FastAPI application."""

from fastapi.testclient import TestClient

from variantagent.api.app import app

client = TestClient(app)


class TestHealthEndpoint:
    def test_health_returns_ok(self) -> None:
        response = client.get("/health")
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "ok"
        assert "version" in data

    def test_health_response_structure(self) -> None:
        response = client.get("/health")
        data = response.json()
        assert set(data.keys()) == {"status", "version"}


class TestAnalyzeEndpoint:
    def test_analyze_returns_not_implemented(self) -> None:
        payload = {
            "sample_id": "test_sample",
        }
        response = client.post("/analyze", json=payload)
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "not_implemented"
