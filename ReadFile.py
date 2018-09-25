"""Proceses a Genbank file and allows for searching, if a pathway exists."""

import sys
import re
from typing import Any, List, Dict, Optional, NamedTuple, Tuple

from collections import namedtuple
from Bio import SeqIO
from CircularGenome import CircularGenome
from Gene import Gene

from custom_blast import blast_seq
import Search as search


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
            cluster_list : List = list()
            pathway_gene : List = self.GENOME.createPathway(parameter, int(search_type), basePairs)
            core_gene = self.GENOME.provide_core_gene()
            blast_output = blast_seq(core_gene)
            self.analyze_other_genomes(blast_output[2], core_gene)
            cluster_list.append(self.get_gene_cluster(pathway_gene))
            return self.get_gene_cluster(pathway_gene)
        else:
            return bool_search

    def get_gene_cluster(self, pathway: List) -> List:
        """Get the necessary info to represent the gene."""
        return [x.serialize() for x in pathway]
    
    def analyze_other_genomes(self, accession_number: str, core_gene: str) -> None:
        """Analyze the results from the blast query."""
        self.GENOME.delete()
        record = search.searchGenbank(accession_number)
        record = SeqIO.read(record, "genbank")
        self.parse_file(record)
        self.GENOME.compare_gene(core_gene)

