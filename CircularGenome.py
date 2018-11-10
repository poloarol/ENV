"""Inserts key-value pairs for genes and their indentifiers and allow for searching, which produced a list of all genes around it by a certain amt of bp."""

import itertools
from collections import namedtuple
from typing import Any, List, Dict, Optional, NamedTuple, Iterator, Tuple

from Bio import Align

import Levenshtein

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
                    self.set_key(main_key)
                    return FOUND
        elif option == self.GENE:
            for main_key in self.genome.keys():
                if value == main_key.gene:
                    self.set_key(main_key)
                    return FOUND
        elif option == self.ID:
            for main_key in self.genome.keys():
                if value == main_key.protein_id:
                    self.set_key(main_key)
                    return FOUND
        else:
            for main_key in self.genome.keys():
                if value == main_key.product:
                    self.set_key(main_key)
                    return FOUND

        return NOT_FOUND

    def createPathway(self, value: str, option: int, basePairs: int) -> List:
        """After haven found a gene exist in the genome, You create a pathway with genes +/- bp."""
        return self.__createPathway(value, option, basePairs)


    def __createPathway(self, value: str, option: int, basePairs: int) -> List:
        """Create the pathway."""
        left = self.__leftPathway(value, option, basePairs)
        right = self.__rightPathway(value, option, basePairs)
        left.reverse()
        left.append(self.genome[self.key])
        left.extend(right)
        return left

    def __rightPathway(self, value: str, option: int, basePairs: int) -> List:
        """Create the rightpathway."""
        pathway: Iterator = itertools.cycle(self.genome)
        if option is self.LOCUS:
            return self.__locus(value, option, basePairs, pathway)
        elif option is self.GENE:
            return self.__gene(value, option, basePairs, pathway)
        elif option is self.ID:
            return self.__protein(value, option, basePairs, pathway)
        else:
            return self.__product(value, option, basePairs, pathway)

    def __leftPathway(self, value: str, option: int, basePairs: int) -> List:
        """Create the leftpathway."""
        get_pathway_keys: List = list(itertools.chain(self.genome))
        # reverse order in order to pass in opposite direction
        get_pathway_keys.reverse()
        reverse_pathway: Iterator = itertools.cycle(get_pathway_keys)
        if option is self.LOCUS:
            return self.__locus(value, option, basePairs, reverse_pathway)
        elif option is self.GENE:
            return self.__gene(value, option, basePairs, reverse_pathway)
        elif option is self.ID:
            return self.__protein(value, option, basePairs, reverse_pathway)
        else:
            return self.__product(value, option, basePairs, reverse_pathway)


    def __locus(self, value: str, option: int, basePairs: int, pathway: Iterator) -> List:
        """Create pathway based on locus."""
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
            elif value == val.locus and started:
                break
            elif value == val.locus:
                started = True
            else:
                pass
        return created_pathway


    def __gene(self, value: str, option: int, basePairs: int, pathway: Iterator) -> List:
        """Create pathway based on gene."""
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
            elif value == val.gene and started:
                break
            elif value == val.gene:
                started = True
            else:
                pass
        return created_pathway

    def __protein(self, value: str, option: int, basePairs: int, pathway: Iterator) -> List:
        """Create pathway based on protein id."""
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
            elif value == val.protein_id and started:
                break
            elif value == val.protein_id:
                started = True
            else:
                pass
        return created_pathway

    def __product(self, value: str, option: int, basePairs: int, pathway: Iterator) -> List:
        """Create pathway based on product name."""
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
            elif value == val.product and started:
                break
            elif value == val.product:
                started = True
            else:
                pass
        return created_pathway

    
    def set_key(self, key: NamedTuple) -> None:
        """Provide easy access to core gene."""
        self.key = key
    
    def delete(self):
        """Empty the dictionary."""
        self.genome.clear()

    def provide_core_gene(self):
        """Provide a fasta format of the core gene."""
        return self.__fastaFormat__()
    
    def __fastaFormat__(self) -> str:
        """Format a seq for blast."""
        return ">{0} \n {1}".format(self.key.gene, self.genome[self.key].translation)

    def compare_gene(self, seq: str, similarity: int) -> List:
        """Build a list of all genes matching the entered query, via Levenshtein string comparison."""
        key_list: List = list()
        SIMILARITY_VALUE = similarity
        for key in self.genome:
            value = round(Levenshtein.ratio(self.genome[key].translation, seq), 2)
            if value >= SIMILARITY_VALUE:
                key_list.append(key)
        # self.gene_separation(key_list, bp)
        return key_list


# def gene_separation(self, my_list: List, bp: int) -> None:
#     """Remove all genes whcih are +/- a certain distance from each other, so as to avoid duplicates."""
#     for i in range(len(my_list)):
#         for j in range(len(my_list)):
#             if my_list[i] == my_list[j] and i != j:
#                 val = abs(self.genome[my_list[i]] - self.genome[my_list[j]])
#                 if val <= bp:
#                     del my_list[j]