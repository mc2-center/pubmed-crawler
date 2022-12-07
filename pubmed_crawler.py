"""PubMed 'Crawler' of MC2 Center Publications.

author: nasim.sanati
maintainer: milen.nikolov
maintainer: verena.chung
"""

import os
import re
import argparse
import getpass
import ssl
from datetime import datetime
import json
import requests

from Bio import Entrez
from bs4 import BeautifulSoup
import synapseclient
import pandas as pd
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows
from openpyxl.styles import Font


def login():
    """Log into Synapse. If env variables not found, prompt user.

    Returns:
        syn: Synapse object
    """
    try:
        syn = synapseclient.login(authToken=os.getenv('SYNAPSE_AUTH_TOKEN'),
                                  silent=True)
    except synapseclient.core.exceptions.SynapseNoCredentialsError:
        print("Credentials not found; please manually provide your",
              "Synapse username and password.")
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
            "interested in only scrapping for new publications."))
    parser.add_argument("-g", "--grant_id",
                        type=str, default="syn21918972",
                        help=("Synapse table/view ID containing grant numbers "
                              "in 'grantNumber' column. (Default: syn21918972)"))
    parser.add_argument("-t", "--table_id",
                        type=str, required=True,
                        help="Current Synapse table holding PubMed info.")
    parser.add_argument("-o", "--output_name",
                        type=str, default="publications_" +
                        datetime.today().strftime('%Y-%m-%d'),
                        help=("Filename for output filename. (Default: "
                              "publications_<current-date>)"))
    return parser.parse_args()


def get_grants(syn, table_id):
    """Get list of grant numbers from dataframe.

    Assumptions:
        Synapse table has `grantNumber`, `consortium`, `theme` columns.

    Returns:
        set: valid grant numbers, e.g. non-empty strings
    """
    print("Querying for grant numbers... ")
    grants = (
        syn.tableQuery(
            f"SELECT grantNumber, consortium, theme FROM {table_id}")
        .asDataFrame()
    )
    print(f"  Number of grants: {len(grants)}\n")
    return grants


def get_pmids(grants):
    """Get list of PubMed IDs using grant numbers as search param.

    Returns:
        set: PubMed IDs
    """
    print("Getting PMIDs from NCBI... ")
    query = " OR ".join(grants['grantNumber'].tolist())
    handle = Entrez.esearch(db="pubmed",
                            term=query,
                            retmax=100_000,
                            retmode="xml",
                            sort="relevance")
    pmids = set(Entrez.read(handle).get('IdList'))
    handle.close()

    # Entrez docs suggests to use HTTP POST when text query is >700
    # characters. If warning is received, replace above code with following:
    # base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    # results = json.loads(requests.post(
    #     f"{base_url}?db=pubmed&term={query}&retmax=100000&retmode=json"))
    # pmids = set(results.get('esearchresult').get('idlist'))

    print(f"  Total unique publications: {len(pmids)}\n")
    return pmids


def parse_grant(grant):
    """Parse for grant number from grant annotation."""
    grant_info = re.search(r"(CA\d+)[ /-]?", grant, re.I)
    return grant_info.group(1).upper()


def get_related_info(pmid):
    """Get related information associated with publication.

    Returns:
        dict: XML results for GEO, SRA, and dbGaP
    """
    handle = Entrez.elink(dbfrom="pubmed",
                          db="gds,sra,gap",
                          id=pmid,
                          retmode="xml")
    results = Entrez.read(handle)[0].get('LinkSetDb')
    handle.close()

    related_info = {}
    for result in results:
        db = re.search(r"pubmed_(.*)", result.get('LinkName')).group(1)
        ids = [link.get('Id') for link in result.get('Link')]
        handle = Entrez.esummary(db=db, id=",".join(ids))
        soup = BeautifulSoup(handle, features="xml")
        handle.close()
        related_info[db] = soup
    return related_info


def parse_geo(info):
    """Parse and return GSE IDs."""
    gse_ids = []
    if info:
        tags = info.find_all('Item', attrs={'Name': "GSE"})
        gse_ids = ["GSE" + tag.text for tag in tags]
    return gse_ids


def parse_sra(info):
    """Parse and return SRX/SRP IDs."""
    srx_ids = srp_ids = []
    if info:
        tags = info.find_all('Item', attrs={'Name': "ExpXml"})
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
        tags = info.find_all('Item', attrs={'Name': "d_study_id"})
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
        'query': query,
        'resultType': "core",
        'format': "json",
        'pageSize': 1_000
    }
    response = json.loads(requests.post(url=pmc_url, data=data).content)
    results = response.get('resultList').get('result')

    grants_list = curr_grants.grantNumber.tolist()
    table = []
    with requests.Session() as session:
        for result in results:
            pmid = result.get('pmid')
            if pmid in pmids:

                # GENERAL INFO
                url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}"
                doi = result.get('doi')
                journal_info = result.get('journalInfo').get('journal')
                journal = journal_info.get(
                    'isoabbreviation', journal_info.get('medlineAbbreviation'))
                year = result.get('pubYear')
                title = result.get('title').rstrip(".")
                authors = [
                    f"{author.get('firstName')} {author.get('lastName')}"
                    for author
                    in result.get('authorList').get('author')
                ]
                abstract = result.get('abstractText', "No abstract available.").replace(
                    "<h4>", " ").replace("</h4>", ": ").lstrip()
                keywords = result.get('keywordList', {}).get('keyword', "")

                # ACCESSIBILITY
                unpaywall_url = f"https://api.unpaywall.org/v2/{doi}?email={email}"
                is_open = json.loads(session.get(unpaywall_url).content).get('is_oa')
                if is_open:
                    accessbility = "Open Access"
                    assay = tissue = tumor_type = ""
                else:
                    accessbility = "Restricted Access"
                    assay = tissue = tumor_type = "Pending Annotation"

                # GRANTS
                grants = result.get('grantList', {}).get('grant', [])
                related_grants = [
                    parse_grant(grant.get('grantId'))
                    for grant in grants
                    if grant.get('grantId')
                    and re.search(r"CA\d", grant.get('grantId'), re.I)
                ]
                related_grants = set(
                    filter(lambda x: x in grants_list, related_grants))

                if related_grants:
                    center = curr_grants.loc[curr_grants['grantNumber'].isin(
                        grants)]
                    consortium = ", ".join(set(center['consortium']))
                    themes = ", ".join(set(center['theme'].sum()))
                else:
                    consortium = themes = ""

                # RELATED INFORMATION
                # Contains: GEO, SRA, dbGaP
                related_info = get_related_info(pmid)
                gse_ids = parse_geo(related_info.get('gds'))
                srx, srp = parse_sra(related_info.get('sra'))
                dbgaps = parse_dbgap(related_info.get('gap'))
                dataset_ids = {*gse_ids, *srx, *srp, *dbgaps}

                # Conslidate all info into a single df, then append to list.
                publication_info = {
                    'Component': ["PublicationView"],
                    'Publication Grant Number': [", ".join(related_grants)],
                    'Publication Consortium Name': [consortium],
                    'Publication Theme Name': [themes],
                    'Publication Doi': [doi],
                    'Publication Journal': [journal],
                    'Pubmed Id': [int(pmid)],
                    'Pubmed Url': [url],
                    'Publication Title': [title],
                    'Publication Year': [int(year)],
                    'Publication Keywords': [", ".join(keywords)],
                    'Publication Authors': [", ".join(authors)],
                    'Publication Abstract': [abstract],
                    'Publication Assay': [assay],
                    'Publication TumorType': [tumor_type],
                    'Publication Tissue': [tissue],
                    'Publication Dataset Alias': [", ".join(dataset_ids)],
                    'Publication Accessibility': [accessbility]
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
        table_name = syn.get(table_id).name
        print(f"Comparing with table: {table_name}...")
        current_pmids = (
            syn.tableQuery(f"SELECT pubMedId FROM {table_id}")
            .asDataFrame()['pubMedId']
            .astype(str)
            .tolist()
        )
        pmids -= set(current_pmids)
        print(f"  New publications found: {len(pmids)}\n")

    if pmids:
        print("Pulling information from publications... ")
        table = pull_info(pmids, grants, email)
    else:
        table = pd.DataFrame()
    print()
    return table


def generate_manifest(syn, table, output):
    """Generate manifest file (xlsx) with given publications data."""
    wb = Workbook()
    ws = wb.active
    ws.title = "manifest"
    for r in dataframe_to_rows(table, index=False, header=True):
        ws.append(r)

    # Get latest CV terms to save as "standard_terms".
    query = ("SELECT attribute, preferredTerm FROM syn26433610 "
             "WHERE attribute <> ''"
             "ORDER BY attribute, preferredTerm")
    cv_terms = syn.tableQuery(query).asDataFrame().fillna("").drop_duplicates()
    ws2 = wb.create_sheet("standard_terms")
    for row in dataframe_to_rows(cv_terms, index=False, header=True):
        ws2.append(row)

    # Style the worksheet.
    ft = Font(bold=True)
    ws2["A1"].font = ft
    ws2["B1"].font = ft
    ws2["C1"].font = ft
    ws2.column_dimensions['A'].width = 18
    ws2.column_dimensions['B'].width = 60
    ws2.column_dimensions['C'].width = 12
    ws2.protection.sheet = True

    wb.save(os.path.join("output", output + ".xlsx"))


def main():
    """Main function."""
    syn = login()
    args = get_args()

    # In order to make >3 Entrez requests/sec, 'email' and 'api_key'
    # params need to be set.
    email = os.getenv('ENTREZ_EMAIL')
    Entrez.email = email
    Entrez.api_key = os.getenv('ENTREZ_API_KEY')

    if not os.environ.get('PYTHONHTTPSVERIFY', '') \
            and getattr(ssl, '_create_unverified_context', None):
        ssl._create_default_https_context = ssl._create_unverified_context

    table = find_publications(syn, args.grant_id, args.table_id.strip(), email)
    if table.empty:
        print("Manifest not generated.")
    else:
        print("Generating manifest... ")

        # Generate manifest with open-access publications listed first.
        generate_manifest(
            syn,
            table.sort_values(by='publicationAccessibility'),
            args.output_name)

    print("-- DONE --")


if __name__ == "__main__":
    main()
