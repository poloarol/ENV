from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation


class Gene():

    """
    Gene:
        Class which stores all the information about a gene from GenBank file
        locus_tag, note, translation, protein_id, product are all set to N/A by
        default and start, stop, codon_start, transl_table and strand are 0
    """

    def __init__(self, locus_tag="N/A", gene="N/A", start=0, stop=0, codon_start=0, table=0, product="N/A", protein_id="N/A", translation="N/A", strand=0):  # noqa
        self.locus_tag = locus_tag
        self.start = start
        self.stop = stop
        self.codon_start = codon_start
        self.table = table
        self.product = product
        self.protein_id = protein_id
        self.translation = translation
        self.strand = strand
        self.gene_name = gene


    def record(self):
        """
            Returns a record of a gene to be written to a genbank file
        """
        sequence = Seq(self.translation, IUPAC.protein)
        rec = SeqRecord(sequence, id=self.protein_id, name=self.locus_tag, \
                        description=self.gene_name)
        loc = SeqFeature(FeatureLocation(self.start, self.stop), \
                         type="CDS", strand=self.strand)
        rec.features.append(loc)
        return rec

    def serialize(self):
        """
            Returns gene information as a json
        """
        return {
            "Gene" : self.gene_name,
            "Product" : self.protein_id,
            "LocusTag" : self.locus_tag,
            "Start" : self.start,
            "Stop" : self.stop,
            "Strand" : self.strand
        }

    def __repr__(self):
        return "<Gene class - Gene Name:{0}, Product:{1}, Locus Tag:{2}> \
        ProteinID: {3} ,Position: ({4}, {5}) Strand:{6}".format(
            self.gene_name, self.product, self.locus_tag, self.protein_id, self.start, self.stop, self.strand
                                                )

    def __str__(self):
        return "Gene: ({0}, {1}, {2}, {3}, {4}, {5})".format(
            self.gene_name, self.product, self.locus_tag, self.start, self.stop, self.strand
                                                )
