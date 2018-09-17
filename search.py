"""Search for genbank file."""

import os
import ssl
import urllib
 
from Bio import Entrez

# add email 
Entrez.email = "adjon081@uottawa.ca"

# with API key, you can make upto 10 request per secs
# We don't have an API key

ssl._create_default_https_context = ssl._create_unverified_context

def searchGenbank(accession_number: str) -> str:
    """Download a genbank record from NCBI."""
    try:
        # record = os.path.expanduser('~/Documents/GitHub/ENV/static/uploads/' + accession_number + ".gb")
        handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="gbwithparts", retmode="txt")
        # local_file = open(record, "w+")
        # local_file.write(handle.read())
        # handle.close()
        # local_file.close()
        return handle
    except urllib.error.HTTPError as error:
        print(error.read())
    return ""
