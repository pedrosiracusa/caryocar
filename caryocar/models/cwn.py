#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Coworking Networks
"""
import networkx
import itertools
from collections import Counter

__author__ = "Pedro Correia de Siracusa"
__copyright__ = "Copyright 2018"

class CWN(networkx.Graph):
    """
    Class for coworking networks. Extends networkx Graph class.
    
    Parameters
    ----------
    cliques : iterable
        An iterable of iterables containing names used to compose cliques 
        in the network.
        
    taxons (optional) : iterable
        An iterable containing taxons recorded by each collector clique. This iterable should be the same length and in the same ordering as the cliques iterable. 
        
    namesMap (optional) : caryocar.NamesMap.
        A caryocar NamesMap object for normalizing nodes names.
        
    Edge attributes
    ---------------
    count: int
        The number of occurrences of the edge.
    
    taxons: list (optional)
        A list with all taxons which were recorded by a pair of collectors (those forming the edge).
    
    weight_hyperbolic: float
        The hyperbolic weight of the edge.
        
    
    Examples
    --------
    
    # Creating a CWN from collectors cliques only
    >>> collectors = [ ['a','b','c'], ['d','e'], ['a','c'] ]
    >>> cwn = CWN(cliques=collectors)
    
    >>> cwn.nodes(data=True)
    { 'a': {'count': 2}, 
      'b': {'count': 1}, 
      'c': {'count': 2}, 
      'd': {'count': 1}, 
      'e': {'count': 1} }    
      
    >>> cwn.edges(data=True)
    [ ('a', 'b', {'count': 1, 'taxons': None, 'weight_hyperbolic': 0.5}), 
      ('a', 'c', {'count': 2, 'taxons': None, 'weight_hyperbolic': 1.5}), 
      ('b', 'c', {'count': 1, 'taxons': None, 'weight_hyperbolic': 0.5}), 
      ('d', 'e', {'count': 1, 'taxons': None, 'weight_hyperbolic': 1.0}) ]
    
    # Creating a CWN with collectors and taxons
    >>> collectors = [ ['a','b','c'], 
                       ['d','e'], 
                       ['a','c'],
                       ['a','c'],
                       ['c','d','e'],
                       ['a','b','c','d'],
                       ['a']]
    >>> taxons = ['t1','t2','t1','t3','t1','t1','t4']    
    >>> cwn = CWN(cliques=collectors,taxons=taxons)
    
    >>> cwn.nodes(data=True)
    { 'b': {'count': 2}, 
      'c': {'count': 5}, 
      'a': {'count': 5}, 
      'd': {'count': 3}, 
      'e': {'count': 2} }
    
    >>> cwn.edges(data=True)
    [ ('b', 'c', {'count': 2, 'taxons': ['t1', 't1'], 'weight_hyperbolic': 0.8333333333333333}), 
      ('b', 'a', {'count': 2, 'taxons': ['t1', 't1'], 'weight_hyperbolic': 0.8333333333333333}), 
      ('b', 'd', {'count': 1, 'taxons': ['t1'], 'weight_hyperbolic': 0.3333333333333333}), 
      ('c', 'a', {'count': 4, 'taxons': ['t1', 't1', 't3', 't1'], 'weight_hyperbolic': 2.8333333333333335}), 
      ('c', 'e', {'count': 1, 'taxons': ['t1'], 'weight_hyperbolic': 0.5}), 
      ('c', 'd', {'count': 2, 'taxons': ['t1', 't1'], 'weight_hyperbolic': 0.8333333333333333}), 
      ('a', 'd', {'count': 1, 'taxons': ['t1'], 'weight_hyperbolic': 0.3333333333333333}), 
      ('d', 'e', {'count': 2, 'taxons': ['t2', 't1'], 'weight_hyperbolic': 1.5})])
    
    """
    def __init__(self, data=None, cliques=None, taxons=None, namesMap=None, **attr):
        """
        Initialization of CWN class.
        
        Parameters
        ----------
        cliques : iterable
            An iterable of iterables containing names used to compose cliques 
            in the network.
            
        namesMap (optional) : caryocar.NamesMap.
            A caryocar NamesMap object for normalizing nodes names.
        """
       
        if cliques is not None:
            if namesMap:
                nmap = namesMap.getMap()
                cliques = [ [ nmap[n] for n in nset ] for nset in cliques ]
            
            # prevent self-loops
            cliques = [ list(set(nset)) for nset in cliques ]
            
            # attributes dicts
            hyperb_weight = lambda ts: 1/(ts-1) 
            e_attr_hyperbWeight=dict()
            e_attr_taxon=dict()
            e_attr_count=dict()
            
            # build edges from records
            if taxons is None: cliques_taxons = map( lambda c: (c,None), cliques)
            else: cliques_taxons = zip(cliques,taxons)
            for clique,taxon in cliques_taxons:
                teamsize=len(clique)
                edgesFromClique = itertools.combinations(clique,2)
                for e in edgesFromClique:
                    e = tuple(sorted(e))
                    e_attr_count[e] = e_attr_count.get(e,0)+1
                    e_attr_taxon[e] = e_attr_taxon.get(e,[])+[taxon] if taxons is not None else None
                    e_attr_hyperbWeight[e] = e_attr_hyperbWeight.get(e,0)+hyperb_weight(teamsize)
            
            edges = e_attr_count.keys()
            data = list(edges)
            
        super().__init__(incoming_graph_data=data,**attr)
    
        # insert nodes and set count attribute
        nodes_counts = Counter( col for clique in cliques for col in clique )
        nodes = nodes_counts.keys()
        
        self.add_nodes_from(nodes)
        networkx.set_node_attributes(self,values=nodes_counts,name='count')
       
        # set edges attributes
        networkx.set_edge_attributes(self,e_attr_count,'count')
        networkx.set_edge_attributes(self,e_attr_taxon,'taxons')
        networkx.set_edge_attributes(self,e_attr_hyperbWeight,'weight_hyperbolic')



if __name__=="__main__":
    
    collectors = [
    # col1, col2, col3 and col4 are connected
    ['col1','col2','col3','col4'],
    ['col1','col2','col3'],
    ['col1','col2','col3'],
    ['col1','col3','col2'],
    ['col1','col2'],
    ['col1','col2'],
    ['col1','col2'],
    ['col1','col3'],
    ['col2','col3'],
    ['col2','col4'],
    ['col2','col4'],
    ['col4'],
    # col5 is isolated
    ['col5'],
    ['col5'],
    # col7 and col8 are connected
    ['col7','col8'],
    ['col7','col8'],
    # col9 would lead to self loop
    ['col9','col9'],
    ['col9','col9'] ]

    cwn = CWN(cliques=collectors)
