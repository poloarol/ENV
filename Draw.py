from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from reportlab.lib.colors import red, grey, orange, green, white
from reportlab.lib.colors import brown, blue, lightblue, purple


class Draw():
    """
        Class used to draw segments of genome
    """
    def __init__(self):
        self.diagram = GenomeDiagram.Diagram("Genome Map")
        self.maxlen = 0

    def draw_gene(self, seq):
        """
            (list, int) -> image (pdf, png, svg)
        """
        col = [red, grey, orange, grey, orange, grey, green, grey, brown, blue, lightblue, grey, purple]  # noqa

        track = self.diagram.new_track(1, name="Annotated Features")
        feature_set = track.new_set()

        i = 0
        color = white

    def split_array(lst):
        """
            (list) -> (list) * len(input)
            Takes a list an splits it into unidimensional arrays
        """
        
