"""Helper method to generate phylogentic trees."""

from Bio import AlignIO, SeqIO, Seq, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import pylab
import os

def alignment(input_file):
    """
    File(.fasta) - > File(.fasta)
    ---------------------------------------------
    This function takes in a DNA or Protein
    and performs a multiple sequence alignment
    and returns a the aligned file
    """

    records = SeqIO.parse(input_file, 'fasta')
    records = list(records) # make a copy, otherwise our generator
                            # is exhausted after calculating maxlen
    maxlen = max(len(record.seq) for record in records)

    # pad sequences so that they all have the same length
    for record in records:
        if len(record.seq) != maxlen:
            sequence = str(record.seq).ljust(maxlen,'.')
            record.seq = Seq.Seq(sequence)
            
    assert all(len(record.seq) == maxlen for record in records)

    # write to temporary file and do all alignment
    output_file = '{}_padded.fasta'.format(os.path.splitext(input_file)[0])
    with open(output_file, 'w') as f:
        SeqIO.write(records, f, 'fasta')
    alignment = AlignIO.read(output_file, 'fasta')

    return make_tree(alignment)

def make_tree(alignment):
    """
    File(.fasta) -> File (.nwk, .nex, .xml)
    ----------------------------------------------------------
    Takes in a fasta file and produces a Tree object to the given file or object
    in a Newick, Nexul or XML form
    Uses upgma method, which works backwards by finding similar organisms and grouping
    them together
    """

    
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignment)

    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)

    tree = Phylo.write(tree, file_name + '.nwk', 'newick')
    return tree

        

# def view_tree(tree_object):
#     """
#     File(.nex,.xml,.nex) -> None
#     ---------------------------------------------------------------------------------------
#     Produces an ascii view of the generated tree
#     """

#     if(tree_object.endswith('.nwk')):
#        t = next(Phylo.parse(tree_object, 'newick'))
#     elif(tree_object.endswith('.xml')):
#        t = next(Phylo.parse(tree_object, 'phyloxml'))
#     else:
#        t = next(Phylo.parse(tree_object, 'nexus'))

    
#     t.rooted = True
#     t.ladderize()

#     Phylo.draw_ascii(t)
#     ##Phylo.draw(t, branch_labels=lambda c: c.branch_length, label_func=get_label)