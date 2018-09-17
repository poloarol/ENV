"""Proceses a Genbank file and allows for searching, if a pathway exists."""

import sys
import re
from typing import Any, List, Dict, Optional, NamedTuple

from collections import namedtuple
from Bio import SeqIO
from CircularGenome import CircularGenome
from Gene import Gene


class ReadFile():
    """Reads a Genbank file and processes it."""

    def __init__(self, filename: str) -> None:
        """Initialize the Circular genome and read the files."""
        self.GENOME: CircularGenome = CircularGenome()
        self.record = SeqIO.read(filename, "genbank")
        self.parse_file(self.record)

    def parse_file(self, record):
        """Parse the file and creates the Genome."""
        locus_tag: str = "N/A"
        product: str = "N/A"
        protein_id: str = "N/A"
        translation: str = "N/A"
        gen: str = "N/A"
        start: int = 0
        stop: int = 0
        codon_start: int = 0
        table: int = 0
        strand: int = 0
        Info: NamedTuple = namedtuple('Info', 'locus, gene, protein_id, product, length')


        for feature in record.features:
            if "CDS" in feature.type:
                for value in feature.qualifiers:
                    if 'protein_id' in feature.qualifiers.keys():
                        protein_id = feature.qualifiers['protein_id'][0]
                    if 'translation' in feature.qualifiers.keys():
                        translation = feature.qualifiers['translation'][0]
                    if 'product' in feature.qualifiers.keys():
                        product = feature.qualifiers['product'][0]
                    if 'locus_tag' in feature.qualifiers.keys():
                        locus_tag = feature.qualifiers['locus_tag'][0]
                    if 'transl_table' in feature.qualifiers.keys():
                        table = feature.qualifiers['transl_table'][0]
                    if 'codon_start' in feature.qualifiers.keys():
                        codon_start = feature.qualifiers['codon_start'][0]
                    if 'gene' in feature.qualifiers.keys():
                        gen = feature.qualifiers['gene'][0]
                strand = feature.strand
                start = feature.location.start.position
                stop = feature.location.end.position
                gene: Gene = Gene(locus_tag, gen, start, stop, codon_start, \
                            table, product, protein_id, translation, strand)
                key: NamedTuple = Info(locus=locus_tag, gene=gen, protein_id=protein_id, \
                           product=product, length=stop-start)
                self.GENOME.add(key, gene)


    def get_gene(self, search_type: str, parameter: str, basePairs: int) -> Any:
        """Search if gene on interest exist at specifed location and create the sub-pathway."""
        bool_search: bool = self.GENOME.findGene(int(search_type), parameter)
        if (bool_search):
            genes: List = self.GENOME.createPathway(parameter, int(search_type), basePairs)
            return [gene.serialize() for gene in genes]
        else:
            return bool_search
