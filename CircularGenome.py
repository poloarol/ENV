from collections import namedtuple

import itertools

class CircularGenome():
    """
        Dictionary used with itertools in order to mimmick the properties of
        a bacterial genome genome. NamedTuples used as keys for the purpose
        of searching a gene
    """

    def __init__(self):
        self.Info = namedtuple('Info', 'locus, gene, protein_id, product, length')
        self.genome = {}

    def add(self, key, value):
        """
            method which add a key and node to the linkedlist
        """
        self.key = self.Info(key.locus, key.gene, key.protein_id, key.product, key.length)
        self.genome[self.key] = value
        self.LOCUS = 1
        self.GENE = 2
        self.ID = 3

    def findGene(self, option, value, basePairs):
        """
            method enables to locate the presence of a gene in the Dictionary
            via the created namd namedtuple
        """

        if option == self.LOCUS:
            for main_key in self.genome.keys():
                print(main_key)
                if value == main_key.locus:
                    return self.__createPathway__(value, main_key, option, basePairs)
        elif option == self.GENE:
            for main_key in self.genome.keys():
                print(main_key)
                if value == main_key.gene:
                    return self.__createPathway__(value, main_key, option, basePairs)
        elif option == self.ID:
            for main_key in self.genome.keys():
                if value == main_key.protein_id:
                    return self.__createPathway__(value, main_key, option, basePairs)
        else:
            for main_key in self.genome.keys():
                if value == main_key.product:
                    return self.__createPathway__(value, main_key, option, basePairs)

        return False

    def __createPathway__(self, key, main_key, option, basePairs):
        left = self.__leftPathway__(key, option, basePairs)
        right = self.__rightPathway__(key, option, basePairs)
        left.append(self.genome[main_key])
        left.extend(right)
        return left

    def __rightPathway__(self, key, option, basePairs):
        pathway = itertools.cycle(self.genome)
        if option is self.LOCUS:
            return self.__locus__(key, option, basePairs, pathway)
        elif option is self.GENE:
            return self.__gene__(key, option, basePairs, pathway)
        elif option is self.ID:
            return self.__protein__(key, option, basePairs, pathway)
        else:
            return self.__product__(key, option, basePairs, pathway)

    def __leftPathway__(self, key, option, basePairs):
        get_pathway_keys = list(itertools.chain(self.genome))
        # reverse order in order to pass in opposite direction
        get_pathway_keys.reverse()
        reverse_pathway = itertools.cycle(get_pathway_keys)
        if option is self.LOCUS:
            return self.__locus__(key, option, basePairs, reverse_pathway)
        elif option is self.GENE:
            return self.__gene__(key, option, basePairs, reverse_pathway)
        elif option is self.ID:
            return self.__protein__(key, option, basePairs, reverse_pathway)
        else:
            return self.__protein__(key, option, basePairs, reverse_pathway)


    def __locus__(self, key, option, basePairs, pathway):
        sum = 0
        started = False
        created_pathway = list()
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


    def __gene__(self, key, option, basePairs, pathway):
        sum = 0
        started = False
        created_pathway = list()
        val = self.Info(None, None, None, None, None)
        while(sum < basePairs):
            val = next(pathway)
            print(val)
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

    def __protein__(self, key, option, basePairs, pathway):
        sum = 0
        started = False
        created_pathway = list()
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

    def __product__(self, key, option, basePairs, pathway):
        sum = 0
        started = False
        created_pathway = list()
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
