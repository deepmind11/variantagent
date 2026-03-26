"""Application configuration via environment variables."""

from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    """VariantAgent configuration.

    All settings can be overridden via environment variables or a .env file.
    """

    # LLM
    llm_provider: str = "anthropic"
    llm_model: str = "claude-sonnet-4-20250514"
    anthropic_api_key: str = ""
    openai_api_key: str = ""

    # NCBI E-utilities
    ncbi_api_key: str = ""
    ncbi_email: str = ""

    # Observability
    langsmith_api_key: str = ""
    langsmith_project: str = "variantagent"
    langsmith_tracing_v2: bool = False

    # Human-in-the-loop
    hitl_confidence_threshold: float = 0.7

    # Logging
    log_level: str = "INFO"

    model_config = {"env_file": ".env", "env_file_encoding": "utf-8"}


settings = Settings()
