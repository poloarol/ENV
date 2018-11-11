"""Proceses a Genbank file and allows for searching, if a pathway exists."""

import sys
import re
from typing import Any, List, Dict, Optional, NamedTuple, Tuple

from collections import namedtuple
from Bio import SeqIO
from CircularGenome import CircularGenome
from Gene import Gene

from custom_blast import blast_seq
import search as search


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
        position: int = 1
        Info: NamedTuple = namedtuple('Info', 'locus, gene, protein_id, product, length')
    
        try:
            organism = record.annotations["source"]
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
                    gene: Gene = Gene(position, organism,locus_tag, gen, start, stop, codon_start, \
                                table, product, protein_id, translation, strand)
                    key: NamedTuple = Info(locus=locus_tag, gene=gen, protein_id=protein_id, \
                            product=product, length=stop-start)
                    self.GENOME.add(key, gene)
                    position += 1
        except:
            # If the file has no CDS field or can't be parsed move to the next
            pass

    def get_gene(self, search_type: str, parameter: str, basePairs: int, similarity: int) -> Any:
        """Search if gene on interest exist at specifed location and create the sub-pathway."""
        bool_search: bool = self.GENOME.findGene(int(search_type), parameter)
        if (bool_search):
            cluster_list : List = list()
            pathway_gene : List = self.GENOME.createPathway(parameter, int(search_type), basePairs)
            cluster_list.append(self.get_gene_cluster(pathway_gene))
            core_gene = self.GENOME.provide_core_gene()
            blast_output = blast_seq(core_gene)
            for i in range(len(blast_output)):
                output: List = self.find_core_gene(blast_output[i], core_gene.split(" ")[2], similarity)
                if(len(output) == 0):
                    pass
                else:
                    try:
                        for i in range(len(output)):
                            self.GENOME.set_key(output[i])
                            ## cluster_list.append(self.get_gene_cluster(pathway_gene))
                            output = self.build_pathway(basePairs)
                            ## output = self.get_gene_cluster(basePairs)
                            cluster_list.append(output)
                            # for x in output:
                            #     cluster_list.append(self.get_gene_cluster(x))
                    except Exception:
                        print(type(output[i]))
            return cluster_list
        else:
            return bool_search

    def get_gene_cluster(self, pathway: List) -> List:
         """Get the necessary info to represent the gene."""
         return [x.serialize() for x in pathway]
    
    def find_core_gene(self, accession_number: str, core_gene: str, similarity: int) -> List:
        """Analyze the results from the blast query and creates a pathway with the found organisms."""
        self.GENOME.delete()
        record = search.searchGenbank(accession_number)
        record = SeqIO.read(record, "genbank")
        self.parse_file(record)
        pathway: List = self.GENOME.compare_gene(core_gene, similarity)
        return pathway
    
    def build_pathway(self, basepairs: int) -> List:
        """Create a new pathway."""
        key = self.GENOME.key
        pathway_gene: List = list()
        try:
            if key.locus is not "N/A":
                value = key.locus
                pathway_gene = self.GENOME.createPathway(value, 1, basepairs)
            elif key.gene is not "N/A":
                value = key.gene
                pathway_gene = self.GENOME.createPathway(value, 2, basepairs)
            elif key.product is not "N/A":
                value = key.protein_id
                pathway_gene = self.GENOME.createPathway(value, 3, basepairs)
            else:
                value = key.product
                pathway_gene = self.GENOME.createPathway(value, 0, basepairs)
            return self.get_gene_cluster(pathway_gene)
        except Exception:
            pass

