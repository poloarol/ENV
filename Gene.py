"""Definition of properties and attibutes of a gene."""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation


class Gene():
    """Stores all necessary information to recognise a gene."""

    def __init__(self,m_id = 0, organism="N/A", locus_tag="N/A", gene="N/A", start=0, stop=0, codon_start=0, table=0, product="N/A", protein_id="N/A", translation="N/A", strand=0):  # noqa
        """Attibutes necessary to make a represent a gene class."""
        self.id = m_id
        self.orgName = organism
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
        """Create a genbank file for a single gene."""
        sequence = Seq(self.translation, IUPAC.protein)
        rec = SeqRecord(sequence, id=self.protein_id, name=self.locus_tag, \
                        description=self.gene_name)
        loc = SeqFeature(FeatureLocation(self.start, self.stop), \
                         type="CDS", strand=self.strand)
        rec.features.append(loc)
        return rec

    def serialize(self):
        """Create a json/dictionary representation of gene."""
        name: str = "";
        if self.locus_tag != "N/A":
            name = self.locus_tag
        elif self.protein_id != "N/A":
            name = self.protein_id
        else:
            name = self.gene_name
        return {
            "id" : self.id,
            "orgName" : self.orgName,
            "start" : self.start,
            "length" : self.stop-self.start,
            "name" : self.gene_name,
            "strand" : "+" if self.strand == 1 else "-",
            "extraclass" : [self.protein_id, self.locus_tag, self.product],
            "AA": self.translation
        }

    def __repr__(self):
        """Class representation of gene for debugging purposes."""
        return "<Gene class - Gene Name:{0}, Product:{1}, Locus Tag:{2}> \
        ProteinID: {3} ,Position: ({4}, {5}) Strand:{6}".format(
            self.gene_name, self.product, self.locus_tag, self.protein_id, self.start, self.stop, self.strand
                                                )

    def __str__(self):
        """Class Human readable representation of string (toString)."""
        return "Gene: ({0}, {1}, {2}, {3}, {4}, {5})".format(
            self.gene_name, self.product, self.locus_tag, self.start, self.stop, self.strand
                                                )
