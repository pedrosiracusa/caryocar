#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Species-collectors Network Module
"""

import networkx
import scipy
import numpy
import copy
from collections import Counter
from sklearn.metrics.pairwise import cosine_similarity

class SpeciesCollectorsNetwork(networkx.Graph):
    """
    Class for Species-collectors networks. Extends networkx Graph class.
    
    Parameters
    ----------
    species : List or iterable
        A list containing names of species to be associated, in order, to elements in the collectors list.
        
    collectors : List or iterable
        A list containing lists of collectors names, to be associated, in order, to elements in the species list.
        
    namesMap (optional) : caryocar.NamesMap
        A caryocar NamesMap object for normalizing nodes names.
    
    Notes
    -----
    For the model to be created both the species and collectors lists must have the same length.
    The ordering of both species and collectors list is important for creating bipartite edges.
    
    Class Methods
    -------------
    .listSpeciesNodes
    .listCollectorsNodes
    .getSpeciesBag
    .getInterestVector
    
    Examples
    --------
    >>> cols=[ ['col1','col2','col3'],
               ['col1','col2'],
               ['col2','col3'],
               ['col4','col5'],
               ['col4'],
               ['col5','col4'] ]
      
    >>> spp=['sp1','sp2','sp3','sp2','sp3','sp2']
    
    >>> scn = SpeciesCollectorsNetwork( species=spp, collectors=cols )
    
    >>> scn.nodes(data=True)
    { 'sp1': {'bipartite': 1, 'count': 1}, 
      'col1': {'bipartite': 0, 'count': 2}, 
      'col2': {'bipartite': 0, 'count': 3}, 
      'col3': {'bipartite': 0, 'count': 2}, 
      'sp2': {'bipartite': 1, 'count': 3}, 
      'sp3': {'bipartite': 1, 'count': 2}, 
      'col4': {'bipartite': 0, 'count': 3}, 
      'col5': {'bipartite': 0, 'count': 2} }
      
    >>> scn.edges(data=True)
    [ ('sp1', 'col1', {'count': 1}), 
      ('sp1', 'col2', {'count': 1}), 
      ('sp1', 'col3', {'count': 1}), 
      ('col1', 'sp2', {'count': 1}), 
      ('col2', 'sp2', {'count': 1}), 
      ('col2', 'sp3', {'count': 1}), 
      ('col3', 'sp3', {'count': 1}), 
      ('sp2', 'col4', {'count': 2}), 
      ('sp2', 'col5', {'count': 2}), 
      ('sp3', 'col4', {'count': 1}) ]    
    """
    def __init__(self, data=None, species=None, collectors=None, namesMap=None, **attr):
        
        self._parseInputData(species,collectors)
        
        self._biadj_matrix = None
        
        set_bipartite_attr=False # a flag for setting bipartite attribute after graph creation
        if species is not None and collectors is not None:
            if namesMap:
                nmap = namesMap.getMap()
                collectors = [ [ nmap[n] for n in nset ] for nset in collectors ]
            
            # build edges
            if len(species)==len(collectors):
                species = list(species)
                collectors = list(collectors)
                
                data = [ (sp,col) for i,sp in enumerate(species) for col in collectors[i] ]
                set_bipartite_attr=True

        super().__init__(data=data,**attr)
        
        if set_bipartite_attr:
            networkx.set_node_attributes( self, values=dict( (n,1) for n in species), name='bipartite' )
            networkx.set_node_attributes( self, values=dict( (n,0) for cols in collectors for n in cols), name='bipartite' )
            
        # set nodes count attribute
        nodes_cols_counts = Counter( c for cols in collectors for c in cols )
        nodes_sp_counts = Counter( species )
        nodes_counts = nodes_cols_counts.copy()
        nodes_counts.update(nodes_sp_counts)
        networkx.set_node_attributes( self, values=nodes_counts, name='count' )

        # set edges count attribute
        edges = data
        edges_counts = Counter(edges)
        for (u,v),cnt in edges_counts.items():
            self[u][v]['count'] = self[u][v].get('count',0)+cnt
            
    
    def _parseInputData( self, species, collectors ):
        # Check format
        if not all( isinstance(lst,list) for lst in collectors ) and \
               all( isinstance(c,str) for lst in collectors for c in lst ):
            raise ValueError("Collectors data input must be in the format of list of lists of strings.")
        
        if not all( isinstance(sp,str) for sp in species ):
            raise ValueError("Species data input must be in the format of list of strings.")
            
        # Check lengths
        if len(species)!=len(collectors):
            raise ValueError("Species and collectors data lists have different lengths.")
        return
    
    def _buildBiadjMatrix( self, col_sp_order=None ):
        if col_sp_order is None:
            col_sp_order=(sorted(self.listCollectorsNodes()),sorted(self.listSpeciesNodes())) 
            
        
        m = networkx.bipartite.biadjacency_matrix(self,
                                                  row_order=col_sp_order[0],
                                                  column_order=col_sp_order[1],
                                                  weight='count')
        
        self._biadj_ix = ( dict( (c,i) for i,c in enumerate(col_sp_order[0]) ), \
                           dict( (s,i) for i,s in enumerate(col_sp_order[1]) ) )
        
        self._biadj_matrix = (*col_sp_order,m)
        
    def _getBiadjMatrix( self ):
        """
        Returns a COPY of the biadjacency matrix
        """
        if self._biadj_matrix is None:
            self._buildBiadjMatrix()
        return copy.deepcopy(self._biadj_matrix)
    
    def remove_nodes_from( self, nodes ):
        """
        Overrides parent method. 
        If nodes removal make isolated nodes those are also removed. 
        """
        super().remove_nodes_from(nodes)
        isolates = list(networkx.isolates(self))
        return super().remove_nodes_from(isolates)
        
    def listSpeciesNodes(self,data=False):
        """
        Lists nodes from the species set.
        
        Parameters
        ----------
        data : string or bool, default=False
            If False only nodes ids are returned.
            If True nodes ids are returned with their respective attribute dicts as (n, attrDict).
            If a string is passed (with an attribute name) then its value is returned in a 2-tuple (n, attrValue).
        
        Returns
        -------
        Either a list of tuples (n,attrDict) or (n,attrValue) where n is the node's id; or a list of nodes id's n.
        
        Note
        ----
        It is not guaranteed that the same order will be mainained in multiple calls of this function.
        """
        spNodes = set( n for n,b in self.nodes(data='bipartite') if b==1 )
        if data==False:
            return [ n for n in self.nodes(data=data) if n in spNodes ]
        else:
            return [ (n,d) for n,d in self.nodes(data=data) if n in spNodes ]
        
    def listCollectorsNodes(self,data=False):
        """
        Lists nodes from the collectors set.
        
        Parameters
        ----------
        data : string or bool, default=False
            If False only nodes ids are returned.
            If True nodes ids are returned with their respective attribute dicts as (n, attrDict).
            If a string is passed (with an attribute name) then its value is returned in a 2-tuple (n, attrValue).
        
        Returns
        -------
        Either a list of tuples (n,attrDict) or (n,attrValue) where n is the node's id; or a list of nodes id's n.
        
        Note
        ----
        It is not guaranteed that the same order will be mainained in multiple calls of this function.
        """
        colNodes = set( n for n,b in self.nodes(data='bipartite') if b==0 )
        if data==False:
            return [ n for n in self.nodes(data=data) if n in colNodes ]
        return [ (n,d) for n,d in self.nodes(data=data) if n in colNodes ]
    
    def getSpeciesBag( self, collector ):
        """
        Parameters
        ----------
        collector : string
          The id of the collector from which to derive the species bag vector.
          
        Returns
        -------
        A tuple (spIds, vector), where the first element is a list containing all species names and
        the second is the vector containing their counts.
        The species bag vector is stored as a 1xn SciPy sparse matrix.
        """
        colList, spList, m = self._getBiadjMatrix()
        i = self._biadj_ix[0][collector]
        vector = m[i]
        return (spList, vector)
    
    def getInterestVector( self, species ):
        """
        Parameters
        ----------
        species : string
          The id of the species from which to derive the interest vector.
          
        Returns
        -------
        A tuple (colIds, vector), where the first element is a list containing all collectors names and
        the second is the vector containing their counts.
        The interest vector is stored as a 1xn SciPy sparse matrix.
        """
        colList, spList, m = self._getBiadjMatrix()
        m = m.transpose()
        i = self._biadj_ix[1][species]
        vector = m[i]
        return (colList,vector)
    
    def _projection_simple_weighting( self, nodesSet, thresh=None ):

        cols,spp,m = self._getBiadjMatrix()
        m.data=numpy.ones(shape=(len(m.data)),dtype=numpy.int)
        g=networkx.Graph()

        if nodesSet=='species': 
            weightsM = scipy.sparse.triu(m.T.dot(m)).tocsr()
            g.add_nodes_from(self.listSpeciesNodes(data=True))
            n=spp

        elif nodesSet=='collectors':
            weightsM = scipy.sparse.triu(m.dot(m.T)).tocsr()
            g.add_nodes_from(self.listCollectorsNodes(data=True))
            n=cols

        else:
            raise ValueError( "nodesSet argument must be 'species' or 'collectors'" )

        weightsM.setdiag(0)
        if thresh is not None:
            weightsM.data = numpy.where( weightsM.data >= thresh, weightsM.data, 0 )
        weightsM.eliminate_zeros()

        for i,row in enumerate(weightsM):
            data=row.data
            colIndices = row.indices
            for j,w in zip(colIndices,data):
                g.add_edge(n[i],n[j],weight=w)

        return g
        
    def _projection_additive_weighting( self,nodesSet, thresh=None ): 
        g = networkx.Graph()
        
        if nodesSet=='species':
                g.add_nodes_from(self.listSpeciesNodes(data=True))
                
        elif nodesSet=='collectors':
                g.add_nodes_from(self.listCollectorsNodes(data=True))
        else:
            raise ValueError( "nodesSet argument must be 'species' or 'collectors'" )
        
        nodes = g.nodes()
        for u in nodes:
            u_nbrs_o = set(self.adj[u]) 
            u_nbrs_i = set( w for v in u_nbrs_o for w in self.adj[v] )-set([u])
            for v in u_nbrs_i:
                v_nbrs_o = set(self.adj[v])
                common_nodes_o = u_nbrs_o & v_nbrs_o
                weight = sum( (self[u][n]['count'] + self[v][n]['count'])/2 for n in common_nodes_o )
                if weight >= thresh:
                    g.add_edge(u,v,weight=weight)        
        
        return g
    
    def _projection_cosine_similarity( self, which, thresh=None ):
        cols,spp,m = self._getBiadjMatrix()
        g = networkx.Graph()
        
        if which=='interest':
            m=m.T
            g.add_nodes_from(self.listSpeciesNodes(data=True))
            n=spp
            
        elif which=='speciesbag':
            g.add_nodes_from(self.listCollectorsNodes(data=True))
            n=cols
            
        else:
            raise ValueError("which argument must be either 'interest' or 'speciesbag'")
            
        simM = scipy.sparse.csr_matrix(cosine_similarity(m))
        simM.setdiag(0)
        if thresh is not None:
            simM.data = numpy.where( simM.data >= thresh, simM.data, 0 )
        simM.eliminate_zeros()
        
        for i,row in enumerate(simM):
            data=row.data
            colIndices=row.indices
            for j,sim in zip(colIndices,data):
                g.add_edge(n[i],n[j],weight=sim)
        
        return g
    
    def project( self, nodesSet, rule='simple_weighting', thresh=None ):
        """
        Generates a SCN projection onto a nodes set, using one of the available rules.
        
        Parameters
        ----------
        nodesSet : str
            The nodes set to project the graph onto. Input can be either 'species' or 'collectors'.
        
        rule : str, default 'simple_weighting'
            The rule that should be used to assign weights to edges in the projected graph. Available rules are: 'simple_weighting', 'additive_weighting', 'cosine_similarity' 
            
        thresh : numerical (optional)
            A weight threshold value for edge creation. If weight value is below threshold the edge is not created.
        """
        if nodesSet not in ['species','collectors']:
            raise ValueError("nodesSet argument must be 'species' or 'collectors'")
            
        if rule=='simple_weighting':
            g=self._projection_simple_weighting(nodesSet=nodesSet,thresh=thresh)
        elif rule=='additive_weighting':
            g=self._projection_additive_weighting(nodesSet=nodesSet,thresh=thresh)
        elif rule=='cosine_similarity':
            which='speciesbag' if nodesSet=='collectors' else 'interest'
            g=self._projection_cosine_similarity(which=which,thresh=thresh)
        else:
            raise ValueError("Invalid projection rule")
            
        return g
        