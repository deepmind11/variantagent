"""FastAPI application for VariantAgent REST API."""

from fastapi import FastAPI
from pydantic import BaseModel

from variantagent import __version__
from variantagent.models.variant import VariantInput

app = FastAPI(
    title="VariantAgent API",
    description="Multi-agent clinical variant interpretation system",
    version=__version__,
)


class HealthResponse(BaseModel):
    status: str
    version: str


@app.get("/health", response_model=HealthResponse)
async def health() -> HealthResponse:
    """Health check endpoint."""
    return HealthResponse(status="ok", version=__version__)


@app.post("/analyze")
async def analyze(variant_input: VariantInput) -> dict:
    """Analyze variants using the multi-agent system.

    Accepts either a VCF file path or manually specified variants.
    Returns a TriageReport with full provenance.
    """
    # TODO: Wire up to orchestrator agent
    return {"status": "not_implemented", "message": "Scaffold only"}
