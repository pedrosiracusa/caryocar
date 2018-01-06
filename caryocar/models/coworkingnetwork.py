#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Coworking Network Module
"""
import networkx
import itertools
from collections import Counter


class CoworkingNetwork(networkx.Graph):
    """
    Class for coworking networks. Extends networkx Graph class.
    
    Parameters
    ----------
    cliques : iterable
        An iterable of iterables containing names used to compose cliques 
        in the network.
        
    namesMap (optional) : caryocar.NamesMap.
        A caryocar NamesMap object for normalizing nodes names.
        
    Class Methods
    -------------
    
    Examples
    --------
    >>> collectors = [ ['a','b','c'], ['d','e'], ['a','c'] ]
    >>> cwn = CoworkingNetwork(cliques=collectors)
    
    >>> cwn.nodes(data=True)
    { 'a': {'count': 2}, 
      'b': {'count': 1}, 
      'c': {'count': 2}, 
      'd': {'count': 1}, 
      'e': {'count': 1} }    
      
    >>> cwn.edges(data=True)
    [ ('a', 'b', {'count': 1}), 
      ('a', 'c', {'count': 2}), 
      ('b', 'c', {'count': 1}), 
      ('d', 'e', {'count': 1}) ]
    
    
    """
    def __init__(self, data=None, cliques=None, namesMap=None, **attr):
        """
        Initialization of CoworkingNetwork class.
        
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
            
            edgesLists = map( lambda n: itertools.combinations(n,r=2), cliques )
            data = [ edge for edgesList in edgesLists for edge in edgesList ]
            
        super().__init__(data=data,**attr)
    
        # insert nodes and set count attribute
        nodes_counts = Counter( col for clique in cliques for col in clique )
        nodes = nodes_counts.keys()
        
        self.add_nodes_from(nodes)
        networkx.set_node_attributes(self,values=nodes_counts,name='count')
        
        # set edges count attribute
        edges = data
        edges_counts = Counter(edges)
        for (u,v),cnt in edges_counts.items():
            self[u][v]['count'] = self[u][v].get('count',0)+cnt