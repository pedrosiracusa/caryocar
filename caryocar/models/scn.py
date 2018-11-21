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

class SCN(networkx.Graph):
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
    
    Methods
    -------------
    .listSpeciesNodes
    .listCollectorsNodes
    .getSpeciesBag
    .getInterestVector
    .remove_nodes_from
    .fromCrsBiadjMatrix
    .project
    .taxonomicAggregation
    
    Examples
    --------
    SCN model initialization
    >>> cols=[ ['col1','col2','col3'],
               ['col1','col2'],
               ['col2','col3'],
               ['col4','col5'],
               ['col4'],
               ['col5','col4'] ]
      
    >>> spp=['sp1','sp2','sp3','sp2','sp3','sp2']
    
    >>> scn = SCN( species=spp, collectors=cols )
    
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
        
        # Class attributes
        self._biadj_matrix = None
        
        # Class construction routine
        if attr.get('initialize_empty')==True:
            super().__init__(data=data,**attr)
            return        
        
        self._parseInputData(species,collectors)

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

        super().__init__(incoming_graph_data=data,**attr)
        
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
    
    @classmethod
    def fromCrsBiadjMatrix( cls, nset1, nset2, m, cols_sp_axes=(0,1) ):
        """
        Creates a SCN network from a scipy CRS biadjacency matrix
        
        Parameters
        ----------
        nset1 : iterable
            Elements that compose nodes set 1 (are the row elements in the crs matrix). Its length must be the same as the number of rows.
        
        nset2 : iterable
            Elements that compose nodes set 2 (are the column elements in the crs matrix). Its length must be the same as the number of columns.
        
        m : scipy CRS sparse matrix
            The biadjacency matrix used to build the SCN.
            
        cols_sp_axes : binary 2-tuple, default (0,1)
            Code to inform whether collectors and species are respectively represented in rows (axis 0) and columns (axis 1) in the biadjacency matrix or the inverse.
            
        Returns
        -------
        A Species Collectors Network
        """
        cols_sp_axes = (1,0) # default (0,1)
        g=cls(initialize_empty=True)

        g.add_nodes_from(nset1, bipartite = 0 if cols_sp_axes==(0,1) else 1)
        g.add_nodes_from(nset2, bipartite = 1 if cols_sp_axes==(0,1) else 0)

        for i,row in enumerate(m):
            data=row.data
            colIndices=row.indices
            for j,cnt in zip(colIndices,data):
                g.add_edge(nset1[i],nset2[j],count=cnt)
        
        return g
    
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
    
    def taxonomicAggregation(self,grouping):
        """
        Generates a taxonomically aggregated version of this network.
        
        Parameters
        ----------
        grouping : dict of iterables (preferrably dict of sets)
            The grouping to be used on aggregation. A valid grouping example is
            >>> grp = { 'family_a': {'sp1','sp2','sp3'},
                        'family_b': {'sp4'},
                        'family_c': {'sp5','sp6'} }
        """
        grouping = dict( (k,set(v)) for k,v in grouping.items() )
        cols,spp,m = self._getBiadjMatrix()
        m=m.T
        ix = self._biadj_ix[1]

        # Clean grouping by removing sp nodes that do not exist in the original network
        itemsToRemove=[]
        spset = set(spp)
        keysToRemove=set()
        for grp,spp in grouping.items():
            itemsToRemove=set()
            for sp in spp:
                if sp not in spset: itemsToRemove.add(sp)
            grouping[grp] -= itemsToRemove
            if len(grouping[grp])==0: keysToRemove.add(grp)  
        for k in keysToRemove: grouping.pop(k)

        # Initialize the new grouping data dictionary
        data = dict()

        # Obtain counts for each new group from species and include in the data dict
        species_counts = dict(self.listSpeciesNodes(data='count'))
        data['count'] = dict( (grp, sum(species_counts[sp] for sp in spp)) for grp,spp in grouping.items() )

        # Create rows to compose the grouping biadjacency matrix
        # Sums up species bags for species in each group
        aggreg_biadj_rows = []
        groups_order = []
        groups_ixes = [ (grp,[ ix[i] for i in ixes ]) for grp,ixes in grouping.items() ]
        for grp,ixes in groups_ixes:
            rows=[]
            for i in ixes:
                rows.append(m.getrow(i))

            grp_matrix = scipy.sparse.vstack(rows,format='csc')
            grp_sum_arr = scipy.sparse.csr_matrix(grp_matrix.sum(axis=0))
            aggreg_biadj_rows.append( (grp,grp_sum_arr) )
            groups_order.append(grp)

        # Create the biadjacency matrix
        aggreg_biadj = ( groups_order, cols, scipy.sparse.vstack([ row for grp,row in aggreg_biadj_rows ], format='csr') )

        # Create the graph from biadj matrix
        g=self.fromCrsBiadjMatrix(*aggreg_biadj,cols_sp_axes=(1,0))

        # set nodes attributes
        networkx.set_node_attributes(g,dict(self.listCollectorsNodes(data='count')),'count')
        networkx.set_node_attributes(g,data['count'],'count')

        return g
    
    def connectedComponentsSubgraphs(self):
        """
        Creates a list of subgraphs with all connected components.
        
        Returns
        -------
        A list of SCN instances, containing each connected component.
        """
        sgs = []
        for sg in networkx.connected_component_subgraphs(self,copy=True):
            scn_i = self.__class__(initialize_empty=True)
            scn_i.add_nodes_from(sg.nodes(data=True))
            scn_i.add_edges_from(sg.edges(data=True))
            sgs.append(scn_i)
            
        return sgs
