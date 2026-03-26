FROM python:3.11-slim

WORKDIR /app

# System dependencies for cyvcf2
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    libc6-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

COPY pyproject.toml .
COPY src/ src/

RUN pip install --no-cache-dir .

COPY data/ data/
COPY streamlit/ streamlit/

EXPOSE 8000 8501

CMD ["uvicorn", "variantagent.api.app:app", "--host", "0.0.0.0", "--port", "8000"]
