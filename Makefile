.PHONY: install dev test lint typecheck format ci clean docker-build docker-up docker-down

install:
	pip install -e .

dev:
	pip install -e ".[dev,streamlit]"
	pre-commit install

test:
	pytest

test-unit:
	pytest tests/unit -v

test-integration:
	pytest tests/integration -v

test-e2e:
	pytest tests/e2e -v

test-cov:
	pytest --cov-report=html && open htmlcov/index.html

lint:
	ruff check src/ tests/

lint-fix:
	ruff check --fix src/ tests/

typecheck:
	mypy src/

format:
	ruff format src/ tests/

ci: lint typecheck test

clean:
	rm -rf dist/ build/ *.egg-info/ .pytest_cache/ .mypy_cache/ htmlcov/ .coverage
	find . -type d -name __pycache__ -exec rm -rf {} +

docker-build:
	docker compose build

docker-up:
	docker compose up -d

docker-down:
	docker compose down

streamlit:
	streamlit run streamlit/app.py

serve:
	uvicorn variantagent.api.app:app --reload --port 8000
