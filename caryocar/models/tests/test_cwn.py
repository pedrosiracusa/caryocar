# -*- coding: utf-8 -*-

import pytest
import networkx
from caryocar.models import CWN
from caryocar.cleaning import NamesMap

@pytest.fixture
def cwn():
    '''A Coworking Network'''
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
    return CWN(cliques=collectors)
   
@pytest.fixture
def cwn_nm():
    '''A Coworking network initialized with a names map'''
    names=[ "col"+str(i) for i in range(1,7) ]
    remapping={ 'col2':'col3',
                'col4':'col5',
                'col3':'COL3',
                'col1':'COL_1' }
    collectors = [ ['col1','col2','col3'],
                   ['col1','col2'],
                   ['col2','col3'],
                   ['col2'],
                   ['col1'],
                   ['col4','col5'],
                   ['col4'] ]
    nm=NamesMap(names=names,normalizationFunc=lambda x: x, remappingIndex=remapping)
    return CWN(cliques=collectors,namesMap=nm)

@pytest.mark.parametrize("col,expectedCount",[
        ('col1',8),
        ('col4',4),
        ('col5',2),
        ('col9',2) ])
def test_cwn_nodes_count_attribute(cwn,col,expectedCount):
    '''Nodes count attribute keeps record of the number of times a collector appears in the dataset'''
    assert(cwn.nodes[col].get('count')==expectedCount)
    
@pytest.mark.parametrize("u,v,expectedCount",[
        ('col1','col2',7),
        ('col2','col4',3),
        ('col7','col8',2)])
def test_cwn_edges_count_attribute(cwn,u,v,expectedCount):
    '''Edges count attribute keeps record of the number of times a tie between collectors appears in the dataset'''
    assert(cwn.edges[(u,v)].get('count')==expectedCount)
    
def test_cwn_isolated_nodes_included(cwn):
    '''Nodes without connections (isolated) are also included in the network'''
    assert 'col5' in cwn.nodes()
    
def test_cwn_selfloops_not_included(cwn):
    '''Self-loops are not included in the network'''
    assert len(cwn['col9'])==0
    
def test_cwn_networkx_subgraphs_works(cwn):
    '''Networkx connected_components_subgraphs function works for superclass'''
    assert all( isinstance(sg,networkx.Graph) for sg in networkx.connected_component_subgraphs(cwn) )

# with names maps
def test_cwn_initialization_namesmap(cwn_nm):
    '''CWN can be initialized with a names map'''
    assert isinstance(cwn_nm,CWN)
    
def test_cwn_remapped_names_exclusion(cwn_nm):
    '''Names which are remapped to others are excluded from the network'''
    assert cwn_nm.nodes().get('col2') is None
    
def test_cwn_remapped_names_addition(cwn_nm):
    '''Counter attribute of names which are target of a remapping get increased'''
    assert cwn_nm.nodes['COL3'].get('count')==4
    
def test_cwn_remapped_names_avoid_selfloops(cwn_nm):
    '''If two names which were in the same record get remapped to each other no self-loop is created'''
    assert len(cwn_nm['col5'])==0
    assert cwn_nm.nodes().get('col4') is None
    
@pytest.mark.parametrize("col,expectedCount",[
        ('COL_1',3),
        ('COL3',4),
        ('col5',2)
        ])
def test_cwn_remapped_names_node_counts(cwn_nm,col,expectedCount):
    '''Nodes count attribute works for remapped names'''
    assert cwn_nm.nodes(data=True)[col].get('count')==expectedCount
    
@pytest.mark.parametrize("u,v,expectedCount",[
        ('COL_1','COL3',2)])
def test_cwn_remapped_names_edge_counts(cwn_nm,u,v,expectedCount):
    '''Edges count attribute works for remapped names'''
    assert cwn_nm.edges[(u,v)].get('count')==expectedCount
    
    
# TODO:
# Test: SCN does not allow 0-degree nodes
# Test: SCN biadj matrix does not change on querying operations
# Test: SCN biadj matrix querying operations retrieve copies of matrix
# Test: CWN hyperbolic weighting in edge attribs.
# Test: CWN taxons in edge attribs.
    
# Execute tests above on script run
if __name__ == '__main__':
    pytest.main(['-v', __file__])