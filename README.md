<h1 align="center">
  MC<sup>2</sup> Center Pubmed Crawler
</h1>

<h3 align="center">
  Publications manifest generator for the 
  <a href="https://cancercomplexity.synapse.org/" target="_blank">Cancer Complexity Knowledge Portal</a> 
  (CCKP)
</h3>
<br/>

<p align="center">
  <img alt="GitHub release (latest by date)" src="https://img.shields.io/github/release/mc2-center/pubmed-crawler?label=latest%20release&display_name=release&style=flat-square">
  <img alt="GitHub Release Date" src="https://img.shields.io/github/release-date/mc2-center/pubmed-crawler?style=flat-square&color=green">
  <img alt="GitHub" src="https://img.shields.io/github/license/mc2-center/pubmed-crawler?style=flat-square&color=orange">
<p>

Manifests for the CCKP can be generated using Docker or Python (3.7+).
Regardless of approach, a One Sage account is required, as well as an
Entrez account (strongly recommended). Failing to provide Entrez credentials
will most likely result in timeout errors from NCBI.

## :whale: Generate with Docker

### Setup

1. Create a file called `.env` and update its contents with your Synapse
   [Personal Access Token] (PAT) and [NCBI account info].

    ```
    # Synapse Credentials
    SYNAPSE_AUTH_TOKEN=<PAT>

    # Entrez Credentials
    ENTREZ_EMAIL=<email>
    ENTREZ_API_KEY=<apikey>
    ```

2. Open a terminal and log into the Synapse Docker registry with your Synapse
    PAT. Once logged in, you should not have to log in again, unless you log
    out or switch Docker registries.

    ```
    docker login docker.synapse.org --username <syn_username>
    ```

    When prompted for a password, enter your PAT.

    You can alternatively log in non-interactively through `STDIN` - this will
    prevent your password from being saved in the shell's history and log files.
    For example, if you saved your PAT into a file called `synapse.token`:

    ```
    cat ~/synapse.token | \
      docker login docker.synapse.org --username <syn_username> --password-stdin
    ```

### Usage

Run the Docker container, replacing `/path/to/.env` with your path to `.env`.

```
docker run --rm -ti \
  --env-file /path/to/.env \
  --volume $PWD/output:/tmp/output:rw \
  docker.synapse.org/syn21498902/pubmed_crawler
```

If this is your first time running the command, Docker will first pull the image
(max. 1-2 minutes) before running the container.

To pull the latest Docker changes, run the following command:

```bash
docker pull docker.synapse.org/syn21498902/pubmed_crawler
```

### Output

Depending on how many new publications have been added to PubMed since the last
scrape (and NCBI’s current requests traffic), this step could take anywhere from
30 seconds to 15ish minutes. Once complete, a manifest will be found in a folder
called `output`, with a name like `publications_manifest_<yyyy-mm-dd>.xlsx`,
where `<yyyy-mm-dd>` is the current date.

## :snake: Generate with Python

### Setup

1. Clone this repo where you want on your local machine, e.g. current directory,
   `Desktop`, etc.

    ```
    git clone https://github.com/mc2-center/pubmed-crawler.git
    ```

2. In the `pubmed-crawler` directory, copy `.envTemplate` as `.env`, then update
   its contents with your Synapse [Personal Access Token] (PAT) and [NCBI account info].

3. Install the dependencies for the Python scripts, ideally in a virtual
   environment, e.g. [conda] or [pyenv]. For example:

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

- [`syn21868591`] is the Synapse table containing publications already curated for the CCKP

PubMed Crawler uses this table to compare against publications found in PubMed,
based on the grant numbers found in the **Portal - Grants Merged** table ([syn21918972]).
To change the table of grants to query PubMed with, use `-g` or `--grantview_id`. For example:

```
python pubmed_crawler.py -t syn21868591 -g syn33657459
```

When using a different table of grants, ensure that its schema has at least the following columns:

- `grantNumber`
- `consortium`
- `theme`

Below is the full usage of the script:

```
usage: pubmed_crawler.py [-h] [-g GRANT_ID] -t TABLE_ID [-o OUTPUT_NAME]

Get PubMed information from a list of grant numbers and put the results into a CSV file.
Table ID can be provided if interested in only scrapping for new publications.

optional arguments:
  -h, --help            show this help message and exit
  -g GRANT_ID, --grant_id GRANT_ID
                        Synapse table/view ID containing grant numbers in 'grantNumber' column. 
                        (Default: syn21918972)
  -t TABLE_ID, --table_id TABLE_ID
                        Current Synapse table holding PubMed info.
  -o OUTPUT_NAME, --output_name OUTPUT_NAME
```

### Output

Any PMIDs found in PubMed that are not found in the Publications table will
be scraped. Depending on the number of new publications (and NCBI’s current
requests traffic), this step could take anywhere from 30 seconds to 15ish
minutes. Once complete, a manifest will be found in a folder called `output`,
with a name like `publications_manifest_<yyyy-mm-dd>.xlsx`, where `<yyyy-mm-dd>`
is the current date.

## :pencil2: Next Steps

Fill out the manifest(s) as needed, using the pre-defined Controlled Vocabulary
listed in **standard_terms** for applicable columns. Once complete, validate
and upload the manifest(s) with the [Data Curator App] (coming soon!).

<!-- Links -->

[synapse account]: https://www.synapse.org/#!RegisterAccount:0
[personal access token]: https://www.synapse.org/#!PersonalAccessTokens:
[ncbi account info]: https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us
[conda]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html
[pyenv]: https://github.com/pyenv/pyenv#getting-pyenv
[data curator app]: https://sagebio.shinyapps.io/csbc_data_curator/
[syn21918972]: https://www.synapse.org/#!Synapse:syn21918972/tables/
[`syn21868591`]: https://www.synapse.org/#!Synapse:syn21868591/tables/
