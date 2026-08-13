"""PubMed 'Crawler' of MC2 Center Publications.

author: nasim.sanati
maintainer: milen.nikolov
maintainer: verena.chung
"""

import argparse
import getpass
import json
import os
import re
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime

import pandas as pd
import requests
import synapseclient
from Bio import Entrez
from bs4 import BeautifulSoup
from http.client import HTTPException
from openpyxl import Workbook
from openpyxl.styles import Font
from openpyxl.utils.dataframe import dataframe_to_rows
from synapseclient.models import Table, query
from urllib.error import HTTPError


def login():
    """Log into Synapse. If env variables not found, prompt user.

    Returns:
        syn: Synapse object
    """
    try:
        syn = synapseclient.login(
            authToken=os.getenv("SYNAPSE_AUTH_TOKEN"), silent=True
        )
    except synapseclient.core.exceptions.SynapseNoCredentialsError:
        print(
            "Credentials not found; please manually provide your",
            "Synapse username and password.",
        )
        username = input("Synapse username: ")
        password = getpass.getpass("Synapse password: ")
        syn = synapseclient.login(username, password, silent=True)
    return syn


def get_args():
    """Set up command-line interface and get arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Get PubMed information from a list of grant numbers and put "
            "the results into a CSV file.  Table ID can be provided if "
            "interested in only scrapping for new publications."
        )
    )
    parser.add_argument(
        "-g",
        "--grant_id",
        type=str,
        default="syn21918972",
        help=(
            "Synapse table/view ID containing grant numbers "
            "in 'grantNumber' column. (Default: syn21918972)"
        ),
    )
    parser.add_argument(
        "-t",
        "--table_id",
        type=str,
        required=True,
        help="Current Synapse table holding PubMed info.",
    )
    parser.add_argument(
        "-o",
        "--output_name",
        type=str,
        default=datetime.today().strftime("%Y-%m-%d") + "_publications-manifest",
        help=(
            "Filename for output filename. (Default: " "publications_<current-date>)"
        ),
    )
    return parser.parse_args()


def get_grants(syn, table_id):
    """Get list of grant numbers from dataframe.

    Assumptions:
        Synapse table has `grantNumber`, `consortium`, `theme` columns.

    Returns:
        set: valid grant numbers, e.g. non-empty strings
    """
    print("Querying for grant numbers... ")
    grants = query(f"SELECT grantNumber, consortium, theme FROM {table_id}")
    print(f"  Number of grants: {len(grants)}\n")
    return grants


def get_pmids(grants):
    """Get list of PubMed IDs using grant numbers as search param.

    Returns:
        set: PubMed IDs
    """
    print("Getting PMIDs from NCBI... ")
    grant_numbers = grants["grantNumber"].tolist()
    query = "[Grant number] OR ".join(grant_numbers) + "[Grant number]"
    handle = Entrez.esearch(
        db="pubmed", term=query, retmax=100_000, retmode="xml", sort="relevance"
    )
    pmids = set(Entrez.read(handle).get("IdList"))
    handle.close()

    # Entrez docs suggests to use HTTP POST when text query is >700
    # characters. If warning is received, replace above code with following:
    # base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    # results = json.loads(requests.post(
    #     f"{base_url}?db=pubmed&term={query}&retmax=100000&retmode=json"))
    # pmids = set(results.get('esearchresult').get('idlist'))

    print(f"  Total unique publications: {len(pmids)}\n")
    return pmids


def parse_grant(pattern, grant):
    """Parse for grant number based on given pattern."""
    grant_info = re.findall(pattern, grant)
    grant_numbers = [
        grant_number.upper().replace(" ", "").replace("/", "").replace("-", "")
        for grant_number in grant_info
    ]
    return grant_numbers


def get_related_info(pmids, batch_size=200, max_retries=3):
    """Get related information for a collection of PMIDs in batched elink calls.

    Network issues may be encountered when making Entrez requests.
    Retry up to `max_retries` times before skipping.

    Returns:
        dict: mapping of pmid -> XML for GEO, SRA, and dbGaP
    """
    result_map = {}
    pmid_list = list(pmids)
    for start in range(0, len(pmid_list), batch_size):
        chunk = pmid_list[start : start + batch_size]
        linksets = []
        for attempt in range(max_retries):
            try:
                handle = Entrez.elink(
                    dbfrom="pubmed",
                    db="gds,sra,bioproject",
                    id=",".join(chunk),
                    retmode="xml",
                )
                linksets = Entrez.read(handle)
                handle.close()
                break
            except (RuntimeError, HTTPException, HTTPError):
                if attempt < max_retries - 1:
                    print(
                        f"  Network issue getting related info for PMID {chunk[0]}..{chunk[-1]}, trying again..."
                    )
                    time.sleep(1)
                else:
                    print(f"  ⚠️ Failed to get related info for PMID {chunk[0]}..{chunk[-1]}. Skipping...")
        for linkset in linksets:
            pmid = str(linkset.get("IdList", [None])[0])
            if pmid is None:
                continue
            related_info = {}
            for link_db in linkset.get("LinkSetDb", []):
                db = re.search(r"pubmed_(.*)", link_db.get("LinkName")).group(1)
                ids = [link.get("Id") for link in link_db.get("Link")]
                handle = Entrez.esummary(db=db, id=",".join(ids))
                soup = BeautifulSoup(handle, features="xml")
                handle.close()
                related_info[db] = soup
            result_map[pmid] = related_info
    return result_map


def _fetch_oa_status(raw_doi, email):
    """Fetch open-access status from Unpaywall for a single DOI.

    Returns:
        tuple: (raw_doi, accessibility)
    """
    if not raw_doi:
        return raw_doi, "Unknown"
    try:
        response = requests.get(
            f"https://api.unpaywall.org/v2/{raw_doi}?email={email}", timeout=10
        )
        response.raise_for_status()
        if response.json().get("is_oa"):
            return raw_doi, "Open Access"
        return raw_doi, "Restricted Access"
    except (requests.exceptions.HTTPError, requests.exceptions.RequestException, json.JSONDecodeError):
        return raw_doi, "Unknown"


def parse_geo(info):
    """Parse and return GSE IDs."""
    gse_ids = []
    if info:
        tags = info.find_all("Item", attrs={"Name": "GSE"})
        gse_ids = ["GSE" + tag.text for tag in tags]
    return gse_ids


def parse_sra(info):
    """Parse and return SRX/SRP IDs."""
    srx_ids = srp_ids = []
    if info:
        tags = info.find_all("Item", attrs={"Name": "ExpXml"})
        srx_ids = [
            re.search(r'Experiment acc="(.*?)"', tag.text).group(1)
            for tag in tags
            if re.search(r'Experiment acc="(.*?)"', tag.text)
        ]
        srp_ids = {
            re.search(r'Study acc="(.*?)"', tag.text).group(1)
            for tag in tags
            if re.search(r'Study acc="(.*?)"', tag.text)
        }
    return srx_ids, srp_ids


def parse_dbgap(info):
    """Parse and return study IDs."""
    gap_ids = []
    if info:
        tags = info.find_all("Item", attrs={"Name": "d_study_id"})
        gap_ids = [tag.text for tag in tags]
    return gap_ids


def pull_info(pmids, curr_grants, email):
    """Create dataframe of publications and their pulled data.

    Returns:
        df: publications data
    """
    pmc_url = "https://www.ebi.ac.uk/europepmc/webservices/rest/searchPOST"
    query = " OR ".join(pmids)
    data = {
        "query": query,
        "resultType": "core",
        "format": "json",
        "pageSize": 1_000,
    }
    response = json.loads(requests.post(url=pmc_url, data=data).content)
    results = response.get("resultList").get("result")

    grants_list = curr_grants.grantNumber.tolist()

    # Filter down to only publications that are in the list of PMIDs and not errata.
    filtered_results = [
        r for r in results
        if r.get("pmid") in pmids
        and "Published Erratum" not in r.get("pubTypeList", {}).get("pubType", [])
    ]
    related_info_map = get_related_info({r.get("pmid") for r in filtered_results})

    # Fetch OA statuses for all qualifying publications.
    unique_dois = {r.get("doi") for r in filtered_results}
    oa_map = {}
    with ThreadPoolExecutor() as executor:
        for raw_doi, accessibility in executor.map(
            lambda doi: _fetch_oa_status(doi, email), unique_dois
        ):
            oa_map[raw_doi] = accessibility

    table = []
    for result in filtered_results:
        pmid = result.get("pmid")

        # GENERAL INFO
        url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}"
        raw_doi = result.get("doi")
        doi = "https://doi.org/" + raw_doi if raw_doi else None
        journal_info = result.get("journalInfo").get("journal")
        journal = journal_info.get(
            "isoabbreviation", journal_info.get("medlineAbbreviation")
        )
        year = result.get("pubYear")
        title = result.get("title").rstrip(".")
        try:
            authors = [
                f"{author.get('firstName')} {author.get('lastName')}"
                for author in result.get("authorList").get("author")
            ]
        except AttributeError:
            authors = []  # There is not an author list with this publication.
        abstract = (
            result.get("abstractText", "No abstract available.")
            .replace("<h4>", " ")
            .replace("</h4>", ": ")
            .lstrip()
        )
        keywords = result.get("keywordList", {}).get("keyword", "")

        # ACCESSIBILITY
        accessbility = oa_map.get(raw_doi, "Unknown")
        if accessbility == "Open Access":
            assay = tissue = tumor_type = ""
        else:
            assay = tissue = tumor_type = "Pending Annotation"

        # GRANTS
        grants = result.get("grantsList", {}).get("grant", [])
        pattern = re.compile(r"(CA[ /-]?\d{6})", re.I)
        related_grants = {
            grant_number
            for grant in grants
            if grant.get("grantId")
            for grant_number in parse_grant(pattern, grant.get("grantId"))
            if grant_number in grants_list
        }
        if related_grants:
            center = curr_grants.loc[
                curr_grants["grantNumber"].isin(related_grants)
            ]
            consortium = ", ".join(set(center["consortium"].sum()))
            themes = ", ".join(set(center["theme"].sum()))
        else:
            consortium = themes = ""

        # RELATED INFORMATION
        # Contains: GEO, SRA, dbGaP
        related_info = related_info_map.get(pmid, {})
        gse_ids = parse_geo(related_info.get("gds"))
        srx, srp = parse_sra(related_info.get("sra"))
        dbgaps = parse_dbgap(related_info.get("gap"))
        dataset_ids = {*gse_ids, *srx, *srp, *dbgaps}

        # Conslidate all info into a single df, then append to list.
        publication_info = {
            "Component": ["PublicationView"],
            "Publication Grant Number": [", ".join(related_grants)],
            "Publication Consortium Name": [consortium],
            "Publication Theme Name": [themes],
            "Publication Doi": [doi],
            "Publication Journal": [journal],
            "Pubmed Id": [int(pmid)],
            "Pubmed Url": [url],
            "Publication Title": [title],
            "Publication Year": [int(year)],
            "Publication Keywords": [", ".join(keywords)],
            "Publication Authors": [", ".join(authors)],
            "Publication Abstract": [abstract],
            "Publication Assay": [assay],
            "Publication Tumor Type": [tumor_type],
            "Publication Tissue": [tissue],
            "Publication Dataset Alias": [", ".join(dataset_ids)],
            "Publication Accessibility": [accessbility],
        }
        row = pd.DataFrame(publication_info)
        table.append(row)
    return pd.concat(table)


def find_publications(syn, grant_id, table_id, email):
    """Get list of publications based on grants of consortia.

    Returns:
        df: publications data
    """
    grants = get_grants(syn, grant_id)
    pmids = get_pmids(grants)

    # If user provided a table ID, only scrape info from publications
    # not already listed in the provided table.
    if table_id:
        table_name = Table(id=table_id).get().name
        id_col = "Pubmed Id" if table_id == "syn52752398" else "pubMedId"
        print(f"Comparing with table: {table_name}...")
        current_pmids = (
            query(f'SELECT "{id_col}" FROM {table_id}')[id_col]
            .astype(str)
            .tolist()
        )
        pmids -= set(current_pmids)
        print(f"  New publications found: {len(pmids)}\n")

    if pmids:
        print("Pulling information from publications... ")
        table = pull_info(pmids, grants, email)
        print(f"  Publications pre-annotated: {len(table.index)}\n")
    else:
        table = pd.DataFrame()
        print()
    return table


def generate_manifest(table, output):
    """Generate manifest file (xlsx) with given publications data."""
    wb = Workbook()
    ws = wb.active
    ws.title = "manifest"
    for r in dataframe_to_rows(table, index=False, header=True):
        ws.append(r)

    # Get latest CV terms to save as "standard_terms".
    annots = ["assay", "tissue", "tumorType"]
    cv_file = "https://raw.githubusercontent.com/mc2-center/data-models/main/all_valid_values.csv"
    cv_terms = pd.read_csv(cv_file)
    cv_terms = cv_terms.loc[
        cv_terms["category"].str.contains("publication")
        | cv_terms["category"].isin(annots)
    ]
    ws2 = wb.create_sheet("standard_terms")
    for row in dataframe_to_rows(cv_terms, index=False, header=True):
        ws2.append(row)

    # Style the worksheet.
    ft = Font(bold=True)
    ws2["A1"].font = ft
    ws2["B1"].font = ft
    ws2.column_dimensions["A"].width = 15
    ws2.column_dimensions["B"].width = 65
    ws2.protection.sheet = True

    wb.save(os.path.join("output", output + ".xlsx"))


def main():
    """Main function."""
    syn = login()
    args = get_args()

    # In order to make >3 Entrez requests/sec, 'email' and 'api_key'
    # params need to be set.
    email = os.getenv("ENTREZ_EMAIL")
    if not email:
        print(
            "⚠️ WARNING: No email address found in the environment.\n"
            "Requests to the Entrez and Unpaywall APIs may be rate-limited or "
            "return incomplete data. For optimal performance, please set the "
            "'ENTREZ_EMAIL' environment variable. See README for more details."
        )

    Entrez.email = email
    Entrez.api_key = os.getenv("ENTREZ_API_KEY")

    table = find_publications(syn, args.grant_id, args.table_id.strip(), email)
    if table.empty:
        print("Manifest not generated.")
    else:
        print("Generating manifest... ")

        # Generate manifest with open-access publications listed first.
        generate_manifest(
            table.sort_values(by="Publication Accessibility"), args.output_name
        )

    print("-- DONE --")


if __name__ == "__main__":
    main()
