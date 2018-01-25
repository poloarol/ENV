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
        self.__locus_tag = locus_tag
        self.__start = start
        self.__stop = stop
        self.__codon_start = codon_start
        self.__table = table
        self.__product = product
        self.__protein_id = protein_id
        self.__translation = translation
        self.__strand = strand
        self.__gene_name = gene

    def set_locus_tag(self, locus):
        """
            (String) -> None
            Assigns a locus tag to the gene
        """
        self.__locus_tag = locus

    def set_codon_start(self, codon_start):
        """
            (String) -> None
            Assigns a codon_start to the gene
        """
        self.__codon_start = codon_start

    def set_trans_table(self, table):
        """
            (String) -> None
            Assigns a translation table to the gene
        """
        self.__table = table

    def set_gene_name(self, gene):
        """
            (String) -> None
            Assigns a gene name to the gene
        """
        self.__gene_name = gene

    def get_length(self):
        """
            (None) -> int
            Returns the length of a gene
        """
        return self.get_stop() - self.get_start()

    def get_locus_tag(self):
        """
            returns the locus tag
        """
        return self.__locus_tag

    def get_start(self):
        """
            returns the start position
        """
        return self.__start

    def get_stop(self):
        """
            returns the stop position
        """
        return self.__stop

    def get_gene_name(self):
        """
            returns the gene name
        """
        return self.__gene_name

    def get_codon_start(self):
        """
            returns the codon_start
        """
        return self.__codon_start

    def get_transl_table(self):
        """
            returns the transl_table
        """
        return self.__table

    def get_product(self):
        """
            returns the product name
        """
        return self.__product

    def get_protein(self):
        """
            returns the protein id
        """
        return self.__protein_id

    def get_translation(self):
        """
            returns the translation
        """
        return self.__translation

    def get_strand(self):
        """
            returns the strand(sense(1) or antisense(-1))
        """
        return self.__strand


    def record(self):
        """
            Returns a record of a gene to be written to a genbank file
        """
        sequence = Seq(self.get_translation(), IUPAC.protein)
        rec = SeqRecord(sequence, id=self.get_protein(), name=self.get_locus_tag(),                     description=self.get_gene_name())  # noqa
        loc = SeqFeature(FeatureLocation(self.get_start(), self.get_stop()), type="CDS", strand=self.get_strand())  # noqa
        rec.features.append(loc)
        return rec

    def serialize(self):
        """
            Returns gene information as a json
        """
        return {
            "Gene" : self.get_gene_name(),
            "Product" : self.get_product(),
            "LocusTag" : self.get_locus_tag(),
            "Start" : self.get_start(),
            "Stop" : self.get_stop(),
            "Strand" : self.get_strand(),
        }
