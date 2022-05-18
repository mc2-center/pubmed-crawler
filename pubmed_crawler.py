"""PubMed Crawler of CSBC/PS-ON Publications.

author: nasim.sanati
author: milen.nikolov
author: verena.chung
"""

import os
import re
import argparse
import getpass
import ssl
from datetime import datetime
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
        syn = synapseclient.login(
            authToken=os.getenv('SYNAPSE_AUTH_TOKEN'),
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
        description="Scrap PubMed information from a list of grant numbers"
        + " and put the results into a CSV file.  Table ID can be provided"
        + " if interested in only scrapping for new publications.")

    # TODO: default to the grants table/view in the "CSBC PS-ON DB" project
    parser.add_argument("-g", "--grantview_id",
                        type=str, default="syn21918972",
                        help="Synapse table/view ID containing grant numbers in"
                        + " 'grantNumber' column. (Default: syn21918972)")
    parser.add_argument("-t", "--table_id",
                        type=str,
                        help="Current Synapse table holding PubMed info.")
    parser.add_argument("-o", "--output_name",
                        type=str, default="publications_"
                        + datetime.today().strftime('%Y-%m-%d'),
                        help="Filename for output filename. (Default:"
                        + " publications_<current-date>)")
    return parser.parse_args()


def convert_to_stringlist(lst):
    """Create StringList of values for Synapse column type.

    Returns:
        str: list of strings
    """
    if lst:
        return "['" + "', '".join(lst) + "']"
    else:
        return "[]"


def make_urls(url, accessions):
    """Create NCBI link for each accession in the iterable.

    Returns:
        str: list of URLs
    """
    url_list = [url + accession for accession in list(accessions)]
    return ", ".join(url_list)


def get_view(syn, table_id):
    """Get Synapse table/data view containing grant numbers.

    Assumptions:
        Syanpse table/view has column called 'grantNumber'

    Returns:
        dataframe: consortiums and their project descriptions.
    """
    results = syn.tableQuery(
        f"select * from {table_id}").asDataFrame()
    return results[~results['grantNumber'].isnull()]


def get_grants(df):
    """Get list of grant numbers from dataframe.

    Assumptions:
        Dataframe has column called 'grantNumber'

    Returns:
        set: valid grant numbers, e.g. non-empty strings
    """
    print("Querying for grant numbers... ", end="")
    grants = set(df.grantNumber.dropna())
    print(f"{len(grants)} found\n")
    return grants


def get_pmids(grants):
    """Get list of PubMed IDs using grant numbers as search param.

    Returns:
        set: PubMed IDs
    """
    print("Getting PMIDs from NCBI...")
    all_pmids = set()

    # Brian's request: add check that pubs. are retreived for each grant number
    count = 1
    for grant in grants:
        handle = Entrez.esearch(db="pubmed", term=grant,
                                retmax=1_000_000, retmode="xml", sort="relevance")
        pmids = Entrez.read(handle).get('IdList')
        handle.close()
        all_pmids.update(pmids)
        count += 1
    print(f"  Total unique publications: {len(all_pmids)}\n")
    return all_pmids


def parse_header(header):
    """Parse header div for pub. title, authors journal, year, and doi."""

    # TITLE
    title = header.find('h1').text.strip()

    # JOURNAL
    journal = header.find('button').text.strip()

    # PUBLICATION YEAR
    pub_date = header.find('span', attrs={'class': "cit"}).text
    year = re.search(r"(\d{4}).*?[\.;]", pub_date).group(1)

    # DOI
    doi_cit = header.find(attrs={'class': "citation-doi"})
    doi = doi_cit.text.strip().lstrip("doi: ").rstrip(".") if doi_cit else ""

    # AUTHORS
    authors = [a.find('a').text for a in header.find_all(
        'span', attrs={'class': "authors-list-item"})]

    return (title, journal, year, doi, authors)


def parse_grant(grant):
    """Parse for grant number from grant annotation."""
    grant_info = re.search(r"(CA\d+)[ /-]", grant, re.I)
    return grant_info.group(1).upper()


def get_related_info(pmid):
    """Get related information associated with publication.

    Entrez will be used for optimal retrieval (since NCBI will kick
    us out if we web-scrap too often).

    Returns:
        dict: XML results for GEO, SRA, and dbGaP
    """
    handle = Entrez.elink(dbfrom="pubmed", db="gds,sra,gap", id=pmid,
                          remode="xml")
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
        tags = info.find_all('item', attrs={'name': "GSE"})
        gse_ids = ["GSE" + tag.text for tag in tags]
    return gse_ids


def parse_sra(info):
    """Parse and return SRX/SRP IDs."""
    srx_ids = srp_ids = []
    if info:
        tags = info.find_all('item', attrs={'name': "ExpXml"})
        srx_ids = [re.search(r'Experiment acc="(.*?)"', tag.text).group(1)
                   for tag in tags]
        srp_ids = {re.search(r'Study acc="(.*?)"', tag.text).group(1)
                   for tag in tags}
    return srx_ids, srp_ids


def parse_dbgap(info):
    """Parse and return study IDs."""
    gap_ids = []
    if info:
        tags = info.find_all('item', attrs={'name': "d_study_id"})
        gap_ids = [tag.text for tag in tags]
    return gap_ids


def scrape_info(pmids, curr_grants, grant_view):
    """Create dataframe of publications and their pulled data.

    Returns:
        df: publications data
    """
    columns = ["doi", "journal", "pubMedId", "bioProjectIds", "bioProjectAccns",
               "pubMedUrl", "publicationTitle", "publicationYear", "keywords",
               "mesh", "authors", "consortium", "grantId", "grantNumber",
               "gseAccns", "gseUrls", "srxAccns", "srxUrls", "srpAccns",
               "srpUrls", "dpgapAccns", "dpgapUrls", "dataset", "tool", "assay",
               "tumorType", "tissue"]

    if not os.environ.get('PYTHONHTTPSVERIFY', '') \
            and getattr(ssl, '_create_unverified_context', None):
        ssl._create_default_https_context = ssl._create_unverified_context

    table = []

    for pmid in pmids:
        session = requests.Session()
        url = f"https://www.ncbi.nlm.nih.gov/pubmed/?term={pmid}"
        soup = BeautifulSoup(session.get(url).content, features="xml")

        if not soup.find(attrs={'aria-label': "500 Error"}):
            # HEADER
            # Contains: title, journal, pub. date, authors, pmid, doi
            header = soup.find(attrs={'id': "full-view-heading"})

            # PubMed utilizes JavaScript now, so content does not always
            # fully load on the first try.
            if not header:
                soup = BeautifulSoup(session.get(url).content, features="xml")
                header = soup.find(attrs={'id': "full-view-heading"})

            title, journal, year, doi, authors = parse_header(header)
            authors = ", ".join(authors)

            # GRANTS
            try:
                grants = [g.text.strip() for g in soup.find(
                    'div', attrs={'id': "grants"}).find_all('a')]

                # Filter out grant annotations not in consortia.
                grants = {parse_grant(grant) for grant in grants
                            if re.search(r"CA\d", grant, re.I)}
                grants = list(filter(lambda x: x in curr_grants, grants))
            except AttributeError:
                grants = []

            # Nasim's note: match and get the grant center Synapse ID from
            # its view table by grant number of this journal study.
            grant_id = consortium = ""
            if grants:
                center = grant_view.loc[grant_view['grantNumber'].isin(
                    grants)]
                grant_id = ", ".join(list(set(center.grantId)))
                consortium = ", ".join(list(set(center.consortium)))

            # KEYWORDS
            abstract = soup.find(attrs={"id": "abstract"})
            try:
                keywords = abstract.find(text=re.compile(
                    "Keywords")).find_parent("p").text.replace(
                        "Keywords:", "").strip()
            except AttributeError:
                keywords = ""

            # MESH TERMS
            mesh = soup.find(attrs={"id": "mesh-terms"})
            try:
                mesh = sorted({term.text.strip().rstrip("*").split(" / ")[0]
                                for term in mesh.find_all(
                    attrs={"class": "keyword-actions-trigger"})})
            except AttributeError:
                mesh = []
            finally:
                mesh = convert_to_stringlist(mesh)

            # RELATED INFORMATION
            # Contains: GEO, SRA, dbGaP
            related_info = get_related_info(pmid)

            gse_ids = parse_geo(related_info.get('gds'))
            gse_url = make_urls(
                "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", gse_ids)

            srx, srp = parse_sra(related_info.get('sra'))
            srx_url = make_urls("https://www.ncbi.nlm.nih.gov/sra/", srx)
            srp_url = make_urls(
                "https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=", srp)

            dbgaps = parse_dbgap(related_info.get('gap'))
            dbgap_url = make_urls(
                "https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=",
                dbgaps)

            row = pd.DataFrame(
                [[doi, journal, int(pmid), "", "", url, title, int(year),
                    keywords, mesh, authors, consortium, grant_id,
                    ", ".join(grants), convert_to_stringlist(
                        gse_ids), gse_url,
                    convert_to_stringlist(srx), srx_url,
                    convert_to_stringlist(list(srp)), srp_url,
                    convert_to_stringlist(dbgaps), dbgap_url,
                    "", "", "", "", ""]],
                columns=columns)
            table.append(row)
        else:
            print(f"{pmid} publication not found - skipping...")
        session.close()

    return pd.concat(table)


def find_publications(syn, grantview_id, table_id):
    """Get list of publications based on grants of consortia.

    Returns:
        df: publications data
    """
    grant_view = get_view(syn, grantview_id)
    grants = get_grants(grant_view)
    pmids = get_pmids(grants)
    table = []

    # If user provided a table ID, only scrape info from publications
    # not already listed in the provided table.
    if table_id:
        print(f"Comparing with table {table_id}...")
        current_publications = syn.tableQuery(
            f"SELECT * FROM {table_id}").asDataFrame()
        current_pmids = {str(pmid)
                         for pmid in list(current_publications.pubMedId)}
        pmids -= current_pmids
        print(f"  New publications found: {len(pmids)}\n")

    if pmids:
        print("Pulling information from publications... ")
        table = scrape_info(pmids, grants, grant_view)
    print("DONE ✓")
    return table


def generate_manifest(syn, table, output):
    """Generate manifest file (xlsx) with given publications data."""
    wb = Workbook()
    ws = wb.active
    ws.title = "manifest"
    for r in dataframe_to_rows(table, index=False, header=True):
        ws.append(r)

    # Get latest CV terms to save as "standard_terms".
    query = ("SELECT key, value, columnType FROM syn26433610 "
             "WHERE key <> '' AND columnType <> '' "
             "ORDER BY key, value")
    cv_terms = syn.tableQuery(
        query).asDataFrame().fillna("").drop_duplicates()
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
    Entrez.email = os.getenv('ENTREZ_EMAIL')
    Entrez.api_key = os.getenv('ENTREZ_API_KEY')

    table = find_publications(syn, args.grantview_id, args.table_id)
    generate_manifest(syn, table, args.output_name)


if __name__ == "__main__":
    main()
