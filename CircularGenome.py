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
        self.key = self.Info(None, None, None, None)

    def add(self, key, value):
        """
            method which add a key and node to the linkedlist
        """
        self.key = self.Info(key.locus, key.gene, key.protein_id, key.product, key.length)
        self.genome[self.key] = value
        self.LOCUS = 1
        self.GENE = 2
        self.ID = 3

    def findGene(self, option, value):
        """
            method enables to locate the presence of a gene in the Dictionary
            via the created namd namedtuple
        """

        if option == self.LOCUS:
            for main_key in self.genome.keys():
                print(main_key)
                if value == main_key.locus:
                    self.key = main_key
                    return True
        elif option == self.GENE:
            for main_key in self.genome.keys():
                print(main_key)
                if value == main_key.gene:
                    self.key = main_key
                    return True
        elif option == self.ID:
            for main_key in self.genome.keys():
                if value == main_key.protein_id:
                    self.key = main_key
                    return True
        else:
            for main_key in self.genome.keys():
                if value == main_key.product:
                    self.key = main_key
                    return True

        return False

    def createPathway(self, value, option, basePairs):
        return self.__createPathway(value, self.key, option, basePairs)


    def __createPathway(self, key, main_key, option, basePairs):
        left = self.__leftPathway(key, option, basePairs)
        right = self.__rightPathway(key, option, basePairs)
        left.append(self.genome[main_key])
        left.extend(right)
        return left

    def __rightPathway(self, key, option, basePairs):
        pathway = itertools.cycle(self.genome)
        if option is self.LOCUS:
            return self.__locus(key, option, basePairs, pathway)
        elif option is self.GENE:
            return self.__gene(key, option, basePairs, pathway)
        elif option is self.ID:
            return self.__protein(key, option, basePairs, pathway)
        else:
            return self.__product(key, option, basePairs, pathway)

    def __leftPathway(self, key, option, basePairs):
        get_pathway_keys = list(itertools.chain(self.genome))
        # reverse order in order to pass in opposite direction
        get_pathway_keys.reverse()
        reverse_pathway = itertools.cycle(get_pathway_keys)
        if option is self.LOCUS:
            return self.__loc(key, option, basePairs, reverse_pathway)
        elif option is self.GENE:
            return self.__ge(key, option, basePairs, reverse_pathway)
        elif option is self.ID:
            return self.__prote(key, option, basePairs, reverse_pathway)
        else:
            return self.__prote(key, option, basePairs, reverse_pathway)


    def __locus(self, key, option, basePairs, pathway):
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


    def __gene(self, key, option, basePairs, pathway):
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

    def __protein(self, key, option, basePairs, pathway):
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

    def __product(self, key, option, basePairs, pathway):
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
