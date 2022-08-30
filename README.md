<h1 align="center">
  MC2 Center Pubmed Crawler
</h1>

<h3 align="center">
  Publications manifest generator for the <a href="https://cancercomplexity.synapse.org/" target="_blank">Cancer Complexity Knowledge Portal</a> (CCKP)
</h3>
<br/>

<p align="center">
  <img src="https://img.shields.io/maintenance/yes/2022?style=flat-square">
  <a href="https://github.com/mc2-center/pubmed-crawler/commits/main"><img src="https://img.shields.io/github/last-commit/mc2-center/pubmed-crawler?color=informational&style=flat-square"></a>
  <a href="https://github.com/mc2-center/pubmed-crawler/issues"><img src="https://img.shields.io/github/issues-raw/mc2-center/pubmed-crawler?color=important&style=flat-square"></a>
<p>

Manifests for the CCKP can be generated using Docker or Python (3.7+).
Regardless of approach, a Synapse account is required, as well as an Entrez
account (strongly recommended).  Failing to provide Entrez credentials will
most likely result in timeout errors from NCBI.

## :whale: Generate with Docker

### Setup

1. Create a file called `.env` and update its contents with your Synapse
[Personal Access Token] (PAT) and [NCBI account info].

```
# Synapse Credentials
SYNAPSE_AUTH_TOKEN=PAT

# Entrez Credentials
ENTREZ_EMAIL=email
ENTREZ_API_KEY=apikey
```

2. Open a terminal and log into the Synapse Docker registry. Provide your
Synapse credentials when prompted.  You should only need to log in once, 
assuming you do not switch registries.

```
docker login docker.synapse.org
```

You can also log in non-interactively through `STDIN` - this will prevent
your password from being saved in the shell's history and log files:

```
cat ~/syn_password.txt | docker login docker.synapse.org --username <syn_username> --password-stdin
```

### Usage

Run the Docker container, replacing `/path/to/.env` with your path to `.env`.

```bash
docker run --rm -ti \
  --env-file /path/to/.env \
  --volume $PWD/output:/tmp/output:rw \
  docker.synapse.org/syn21498902/pubmed_crawler
```

If this is your first time running the command, Docker will first pull the image
(max. 1-2 minutes) before running the container.

Depending on how many new publications have been added to PubMed since the last
scrape (and NCBI’s current requests traffic), this step could take anywhere from
30 seconds to 15ish minutes.  Once complete, a manifest will be found in a folder
called `output`, with a name like `publications_manifest_<yyyy-mm-dd>.xlsx`, as
well as manifest templates for datasets, files, and tools (where `<yyyy-mm-dd>`
is the current date).

## :snake: Generate with Python

### Setup

1. Clone this repo where you want on your local machine, e.g. current directory,
`Desktop`, etc.

```bash
git clone https://github.com/mc2-center/pubmed-crawler.git
```

2. In the `pubmed-crawler` directory, copy `.envTemplate` as `.env`, then update
its contents with your Synapse [Personal Access Token] (PAT) and [NCBI account info].

3. Install the dependencies for the Python scripts, ideally in a virtual
environment, e.g. [conda] or [pyenv].  For example:
  
```
conda create -n pubmed-crawler python=3.9
conda activate pubmed-crawler
pip install -r requirements.txt
```
  
4. Set environment variables from `.env` so that the scripts will have access
to the credentials.

```
export $(grep -v '^#' .env | xargs)
```

### Usage

While in the virtual environment, run the command:

```
python pubmed_crawler.py -t syn21868591
```

where:

* [`syn21868591`] is the Synapse table containing publications already curated for the CCKP

PubMed Crawler uses this table to compare against publications found in PubMed, 
based on the grant numbers found in the **Portal - Grants Merged** table ([syn21918972]).
To change the table of grants to compare against, use `-g` or `--grantview_id`:

```
python pubmed_crawler.py -t syn21868591 -g syn123
```

Any PMIDs found in PubMed that are not found in the Publications table will
be added to the output file, `publications_manifest_<yyyy-mm-dd>.xlsx`, in the
`output` folder, where `<yyyy-mm-dd>` is the current date.

If there are accompanying datasets, data files, and/or tools, be sure to generate
their templates as well (templates will be found in `output`):

```
python generate_templates.py
```

## :pencil2: Next Steps
Fill out the manifest(s) as needed, using the pre-defined Controlled Vocabulary
listed in **standard_terms** for applicable columns.  Once complete, validate
and upload the manifest(s) with the [MC2 Data Curator App].



<!-- Links -->
[Synapse account]: https://www.synapse.org/#!RegisterAccount:0
[Personal Access Token]: https://www.synapse.org/#!PersonalAccessTokens:
[NCBI account info]: https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us
[conda]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html
[pyenv]: https://github.com/pyenv/pyenv#getting-pyenv
[MC2 Data Curator App]: https://sagebio.shinyapps.io/csbc_data_curator/
[syn21918972]: https://www.synapse.org/#!Synapse:syn21918972/tables/
[`syn21868591`]: https://www.synapse.org/#!Synapse:syn21868591/tables/
