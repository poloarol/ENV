"""Search for genbank file."""

import os
import ssl
from Bio import SeqIO, Entrez

Entrez.email = "adjon081@uottawa.ca"
ssl._create_default_https_context = ssl._create_unverified_context

def searchGenbank(accession_number: str) -> None:
    """Download a genbank record from NCBI."""
    record = accession_number + ".gb"
    handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="gbwithparts", retmode="txt")
    local_file = open(record, "w+")
    local_file.write(handle.read())
    handle.close()
    local_file.close()

def main():
    """Call search method."""
    searchGenbank("AE017125")


if __name__ == "__main__":
    main()