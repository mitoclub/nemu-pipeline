import requests
import time
from xml.etree import ElementTree

def get_pubmed_titles(query, max_results=10):
    # Define the base URL for ESearch
    base_esearch_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'

    # Define the base URL for EFetch
    base_efetch_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'

    # Set the parameters for the ESearch
    esearch_params = {
        'db': 'pubmed',      # Database to search
        'term': query,       # Search term/phrase
        'retmode': 'xml',    # Return mode
        'retmax': max_results  # Number of results to return
    }

    # Perform the ESearch to get the list of PMIDs
    esearch_response = requests.get(base_esearch_url, params=esearch_params)
    esearch_tree = ElementTree.fromstring(esearch_response.content)

    # Extract PMIDs from the ESearch result
    id_list = [id_tag.text for id_tag in esearch_tree.findall('./IdList/Id')]

    # If no IDs were found, return empty list
    if not id_list:
        return []

    # Set the parameters for the EFetch
    efetch_params = {
        'db': 'pubmed',      # Database to fetch from
        'id': ','.join(id_list),  # PMIDs to fetch
        'retmode': 'xml',    # Return mode
        'rettype': 'abstract'  # Return type
    }

    # Perform the EFetch to get the article details
    efetch_response = requests.get(base_efetch_url, params=efetch_params)
    print(str(efetch_response.content).split('ArticleTitle')[1])

with open('ref_data.txt') as fin:
    i = 0
    for line in fin:
        if line == '\n':
            continue
        i += 1
        doi_link, pubmed_link = line.strip().split()
        doi = doi_link.strip().split('https://doi.org/')[1]

        print(i)
        print(doi_link)
        print(pubmed_link)
        get_pubmed_titles(doi)
        print()
        time.sleep(1)
