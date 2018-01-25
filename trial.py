from collections import namedtuple

from Bio import SeqIO

from CircularGenome import CircularGenome
from Gene import Gene

def paran(self, my_list):
    ret = ""
    for i in range(len(my_list)):
        if my_list[i] != '[' or my_list[i] != ']':
            ret += my_list[i]
    return ret

def parse_file(filename):
    record = SeqIO.read(filename, "genbank")
    locus_tag = product = protein_id = translation = gen = "N/A"
    start = stop = codon_start = table = strand = 0
    Info = namedtuple('Info', 'locus, gene, protein_id, product, length')


    for feature in record.features:
        if "CDS" in feature.type:
            for value in feature.qualifiers:
                if 'protein_id' in feature.qualifiers.keys():
                    protein_id = self.paran(feature.qualifiers['protein_id'])
                if 'translation' in feature.qualifiers.keys():
                    translation = self.paran(feature.qualifiers['translation'])
                if 'product' in feature.qualifiers.keys():
                    product = self.paran(feature.qualifiers['product'])
                if 'locus_tag' in feature.qualifiers.keys():
                    locus_tag = self.paran(feature.qualifiers['locus_tag'])
                if 'transl_table' in feature.qualifiers.keys():
                    table = self.paran(feature.qualifiers['transl_table'])
                if 'codon_start' in feature.qualifiers.keys():
                    codon_start = self.paran(feature.qualifiers['codon_start'])
                if 'gene' in feature.qualifiers.keys():
                    gen = self.paran(feature.qualifiers['gene'])
            strand = feature.strand
            start = feature.location.start.position
            stop = feature.location.end.position
            gene = Gene(locus_tag, gen, start, stop, codon_start, table,product, protein_id, translation, strand)  # noqa
            data = Info(locus=locus_tag, gene=gen, protein_id=protein_id,product=product, length=stop-start)  # noqa
            GENOME.add_node(data, gene)  # noqa


def get_gene(self, search_type, parameter, bases):
    gene_set = GENOME.giveNode(search_type, parameter, bases)
    c = []
    for i in range(len(gene_set)):
        c.append(gene_set[i].serialize())
        print(gene_set[i].serialize(), sys.stdout)
    return c

parse_file()
