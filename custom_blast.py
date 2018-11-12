"""Perform blast on a given gene."""
import ssl

from typing import List

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

ssl._create_default_https_context = ssl._create_unverified_context

def blast_seq(seq: str, expect: int = 10, hitlist_size: int = 100) -> List:
    """Perform a blast on a given gene."""
    result_handle = NCBIWWW.qblast("tblastn", "nr", seq, format_type="XML", hitlist_size=hitlist_size, expect=expect, service="psi")
    blast_records = NCBIXML.parse(result_handle)
    access_num_list = list()
    for rec in blast_records:
        for alignment in rec.alignments:
            access_num_list.append(alignment.accession)
    access_num_list.reverse()
    return access_num_list