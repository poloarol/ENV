"""Inserts key-value pairs for genes and their indentifiers and allow for searching, which produced a list of all genes around it by a certain amt of bp."""

import itertools
from collections import namedtuple
from typing import Any, List, Dict, Optional, NamedTuple, Iterator, Tuple

from Bio import Align

class CircularGenome():
    """
    Dictionary used with itertools in order to mimmick the properties of a bacterial genome genome.
    
    NamedTuples used as keys for the purpose of searching a gene.
    """

    def __init__(self):
        """Initialize the genome with its key values when called."""
        self.Info = namedtuple('Info', 'locus, gene, protein_id, product, length')
        self.genome: Dict = {}
        self.key: Info = self.Info(None, None, None, None, None)
        self.LOCUS: int = 1
        self.GENE: int = 2
        self.ID: int = 3

    def add(self, key, value) -> None:
        """Add a new item to the dictionary."""
        cle = self.Info(key.locus, key.gene, key.protein_id, key.product, key.length)
        self.genome[cle] = value
    

    def findGene(self, option: int, value: str) -> bool:
        """Enable to locate the presence of a gene in the Dictionary via the created namd namedtuple."""
        FOUND = True
        NOT_FOUND = False
        if option == self.LOCUS:
            for main_key in self.genome.keys():
                if value == main_key.locus:
                    self.key = main_key
                    return FOUND
        elif option == self.GENE:
            for main_key in self.genome.keys():
                if value == main_key.gene:
                    self.key = main_key
                    return FOUND
        elif option == self.ID:
            for main_key in self.genome.keys():
                if value == main_key.protein_id:
                    return FOUND
        else:
            for main_key in self.genome.keys():
                if value == main_key.product:
                    self.key = main_key
                    return FOUND

        return NOT_FOUND

    def createPathway(self, value: str, option: int, basePairs: int) -> List:
        """After haven found a gene exist in the genome, You create a pathway with genes +/- bp."""
        return self.__createPathway(value, option, basePairs)


    def __createPathway(self, key: str, option: int, basePairs: int) -> List:
        """Create the pathway."""
        left = self.__leftPathway(key, option, basePairs)
        right = self.__rightPathway(key, option, basePairs)
        left.append(self.genome[self.key])
        left.extend(right)
        return left

    def __rightPathway(self, key: str, option: int, basePairs: int) -> List:
        """Create the rightpathway."""
        pathway: Iterator = itertools.cycle(self.genome)
        if option is self.LOCUS:
            return self.__locus(key, option, basePairs, pathway)
        elif option is self.GENE:
            return self.__gene(key, option, basePairs, pathway)
        elif option is self.ID:
            return self.__protein(key, option, basePairs, pathway)
        else:
            return self.__product(key, option, basePairs, pathway)

    def __leftPathway(self, key: str, option: int, basePairs: int) -> List:
        """Create the leftpathway."""
        get_pathway_keys: List = list(itertools.chain(self.genome))
        # reverse order in order to pass in opposite direction
        get_pathway_keys.reverse()
        reverse_pathway: Iterator = itertools.cycle(get_pathway_keys)
        if option is self.LOCUS:
            return self.__locus(key, option, basePairs, reverse_pathway)
        elif option is self.GENE:
            return self.__gene(key, option, basePairs, reverse_pathway)
        elif option is self.ID:
            return self.__protein(key, option, basePairs, reverse_pathway)
        else:
            return self.__product(key, option, basePairs, reverse_pathway)


    def __locus(self, key: str, option: int, basePairs: int, pathway: Iterator) -> List:
        sum: int = 0
        started: bool = False
        created_pathway: List = list()
        val = self.Info(None, None, None, None, None)
        while(sum < basePairs):
            val = next(pathway)
            if started:
                created_pathway.append(self.genome[val])
                sum += val.length
            elif sum >= basePairs:
                break
            elif key == val.locus and started:
                break
            elif key == val.locus:
                started = True
            else:
                pass
        return created_pathway


    def __gene(self, key: str, option: int, basePairs: int, pathway: Iterator) -> List:
        sum: int = 0
        started: bool = False
        created_pathway: List = list()
        val = self.Info(None, None, None, None, None)
        while(sum < basePairs):
            val = next(pathway)
            if started:
                created_pathway.append(self.genome[val])
                sum += val.length
            elif sum >= basePairs:
                break
            elif key == val.gene and started:
                break
            elif key == val.gene:
                started = True
            else:
                pass
        return created_pathway

    def __protein(self, key: str, option: int, basePairs: int, pathway: Iterator) -> List:
        sum: int = 0
        started: bool = False
        created_pathway: List = list()
        val = self.Info(None, None, None, None, None)
        while(sum < basePairs):
            val = next(pathway)
            if started:
                created_pathway.append(self.genome[val])
                sum += val.length
            elif sum >= basePairs:
                break
            elif key == val.protein_id and started:
                break
            elif key == val.protein_id:
                started = True
            else:
                pass
        return created_pathway

    def __product(self, key: str, option: int, basePairs: int, pathway: Iterator) -> List:
        sum: int = 0
        started: bool = False
        created_pathway: List = list()
        val = self.Info(None, None, None, None, None)
        while(sum < basePairs):
            val = next(pathway)
            sum += val.length
            if started:
                created_pathway.append(self.genome[val])
                sum += val.length
            elif sum >= basePairs:
                break
            elif key == val.product and started:
                break
            elif key == val.product:
                started = True
            else:
                pass
        return created_pathway
    
    def delete(self):
        """Empty the dictionary."""
        self.genome.clear()

    def provide_core_gene(self):
        """Provide a fasta format of the core gene."""
        return self.__fastaFormat__()
    
    def __fastaFormat__(self) -> str:
        """Format a seq for blast."""
        return ">{0} \n{1}".format(self.key.gene, self.genome[self.key].translation)

    def compare_gene(self, seq: str) -> None:
        """Build a list of all genes matching the entered query, via pairwise alignment."""
        alignment = "alignment.txt"
        aligner = Align.PairwiseAligner()
        with open(alignment, "w+") as f:
            for key in self.genome:
                align = aligner.align(self.genome[key].translation, seq)
                for x in sorted(align):
                    f.write(x, x.score)
        f.close()
