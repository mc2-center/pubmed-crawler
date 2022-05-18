<h1 align="center">
  CSBC/PS-ON Pubmed Crawler
</h1>

<h3 align="center">
  Publications manifest generator for the <a href="https://staging.cancercomplexity.synapse.org/" target="_blank">Cancer Complexity Knowledge Portal</a> (CCKP)
</h3>
<br/>

<p align="center">
  <img src="https://img.shields.io/maintenance/yes/2022?style=flat-square">
  <a href="https://github.com/mc2-center/pubmed-crawler/commits/main"><img src="https://img.shields.io/github/last-commit/mc2-center/pubmed-crawler?color=informational&style=flat-square"></a>
  <a href="https://github.com/mc2-center/pubmed-crawler/issues"><img src="https://img.shields.io/github/issues-raw/mc2-center/pubmed-crawler?color=important&style=flat-square"></a>
<p>

Manifests for the CCKP can be generated using [Docker] or Python (3.7+).
Regardless of approach, a Synapse account is required, as well as an Entrez
account (strongly recommended).  Failing to provide Entrez credentials will
most likely result in timeout errors from NCBI.

## Generate with Docker

### Setup

1. Create a file called `.env` and update its contents with your Synapse
[Personal Access Token] (PAT) and [NCBI account info].

```
# Synapse Credentials
SYNAPSE_AUTH_TOKEN="PAT"

# Entrez Credentials
ENTREZ_EMAIL=email
ENTREZ_API_KEY=apikey
```

2. Open a terminal and log into the Synapse Docker hub.

```
docker login docker.synapse.org
```

### Usage

Run the Docker container, replacing `/path/to/.env` with your path to `.env`.

```bash
docker run --rm -ti \
  --env-file /path/to/.env \
  --volume $PWD/output:/tmp/output:rw \
  docker.synapse.org/syn7080714/pubmed_crawler:v2.1.0
```

If this is your first time running the command, Docker will first pull the image
(max. 1-2 minutes) before running the container.

Depending on how many new publications have been added to PubMed since the last
scrape (and NCBI’s current requests traffic), this step could take anywhere from
10 seconds to 15ish minutes.  Once complete, a manifest will be found in a folder
called `output`, with a name like `publications_manifest_yyyy-mm-dd.xlsx`, as
well as manifest templates for datasets, files, and tools.

## Generate with Python

### Setup

1. Download `pubmed-crawler` where you want on your local machine, e.g. `$HOME`,
`Desktop`, etc.

```bash
git clone https://github.com/mc2-center/pubmed-crawler.git
```

2. In the `pubmed-crawler` directory, copy `.envTemplate` as `.env`, then update
its contents with your Synapse [Personal Access Token] (PAT) and [NCBI account info].

3. Install the dependencies for the Python scripts, ideally in a virtual
environment, e.g. [conda] or [pyenv].

```
pip install -r requirements.txt
```

4. Set environment variables from `.env` so that the scripts will have access
to the credentials.

```
export $(grep -v '^#' .env | xargs)
```

### Usage

Run the command:

```
python pubmed_crawler.py -t syn21868591
```

This will pull all grant numbers from **Portal - Grants Merged** (syn21918972),
use them as the search terms in PubMed, then compare the PMIDs found in PubMed
with the existing ones in **Portal - Publications Merged** (syn21868591). The
final output will be a new manifest named `publications_manifest_yyyy-mm-dd.xlsx`
in the `output` folder.

If there are accompanying datasets, data files, and/or tools, be sure to generate
their templates as well:

```
python generate_templates.py
```

## Next Steps
Fill out the manifest(s) as needed, using the pre-defined Controlled Vocabulary
listed in **standard_terms** for applicable columns.  Once complete, validate
upload the manifest(s) with the MC2 Data Curator App (coming soon!).

[Docker]: https://www.docker.com/get-started
[Personal Access Token]: https://www.synapse.org/#!PersonalAccessTokens:
[NCBI account info]: https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us
[conda]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html
[pyenv]: https://github.com/pyenv/pyenv#getting-pyenv