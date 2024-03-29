FROM python:3.9.1-slim-buster

WORKDIR /tmp

COPY requirements.txt .
RUN pip install --upgrade pip \
    && pip install -r requirements.txt
COPY . .

ENTRYPOINT [ "python", "pubmed_crawler.py", "-t syn21868591" ]