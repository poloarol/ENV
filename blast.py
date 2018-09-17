
import ssl

from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIXML

ssl._create_default_https_context = ssl._create_unverified_context

record = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF", 
                    IUPAC.protein),
                    id="YP_025292.1", name="HokC",
                    description="toxic membrane protein")

result_handle = NCBIWWW.qblast("blastp", "nt", record.seq)