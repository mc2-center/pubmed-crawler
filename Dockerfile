FROM python:3.12-slim

COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

WORKDIR /tmp

COPY pyproject.toml .
RUN uv sync --no-dev --no-install-project

COPY . .

ENTRYPOINT ["uv", "run", "pubmed_crawler.py", "-t", "syn21868591"]
