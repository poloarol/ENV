from collections import namedtuple


class CircularGenome():
    """
        Circular bidirectional linkedlist to simulate bacterial genome
    """

    class __Node(object):
        """
            Node class
        """

        __obj = 'details', '_data', '_prev', '_next'
        __info = '_locus', '_gene', '_product', '_protein_id', '_len'

        def __init__(self, Info, data, prev, next):
            self._locus = Info.locus
            self._gene = Info.gene
            self._protein_id = Info.protein_id
            self._product = Info.product
            self._len = Info.length
            self._data = data
            self._prev = prev
            self._next = next

    def __init__(self):
        Info = namedtuple('Info', 'locus, gene, protein_id, product, length')
        # Info serves a key
        self._head = self.__Node(Info, None, None, None)
        self._tail = self.__Node(Info, None, None, None)
        self._head._next = self._tail
        self._tail._prev = self._head
        self._size = 0

    def _behind(self):
        """
            returns the previous element in the list. If the previous is
            the head, returns the head.previous
        """
        if(self._prev is self._head):
            return self._head._prev
        return self._prev

    def add_node(self, d, e):
        """
            method which add a key and node to the linkedlist
        """
        global p
        if self.is_empty():
            p = self.__Node(d, e, self._head, self._head)
            self._head._next = p
            self._head._prev = p
        else:
            p = self._head
            while p._next is not self._head:
                p = p._next
            p._next = self.__Node(d, e, p, self._head)
            self._head._prev = p._next
        self._size += 1
        return p._data

    def is_empty(self):
        """
         checks if the list is empty
        """
        return self.size() == 0

    def size(self):
        """
         gives the size of the list
        """
        return self._size

    def giveNode(self, type, s, b):
        """
            returns a list of all neighbouring genes within specified range
        """
        node_list = []

        if type is 1:
            p = self.__search_by_locus(s)
        elif type is 2:
            p = self.__search_by_gene(s)
        elif type is 3:
            p = self.__search_by_id(s)
        else:
            p = self.__search_by_prod(s)

        if p is False:
            return "Item can't be found"

        p_ahead = p_behind = p
        length_ahead = length_behind = 0

        while length_ahead <= b:
            p_ahead = p_ahead._next
            if p_ahead is self._head:
                p_ahead = self._head._next
            node_list.append(p_ahead._data)
            length_ahead += p_ahead._len

        while length_behind <= b:
            p_behind = p_behind._prev
            if p_behind is self._head:
                p_behind = self._head._prev
            node_list.append(p_behind._data)
            length_behind += p_behind._len

        node_list.append(p._data)

        return node_list

    def __search_by_locus(self, s):
        p = self._head._next
        while p._locus != s:
            p = p._next
            if p is self._head:
                return False
        return p

    def __search_by_gene(self, s):
        p = self._head._next
        while p._gene != s:
            p = p._next
            if p is self._head:
                return False
        return p

    def __search_by_id(self, s):
        p = self._head._next
        while p._protein != s:
            p = p._next
            if p is self._head:
                return False
        return p

    def __search_by_prod(self, s):
        p = self._head._next
        while p._product != s:
            p = p._next
            if p is self._head:
                return False
        return p
