# -*- coding: utf-8 -*-

import pytest
import networkx
from caryocar.models import CoworkingNetwork


@pytest.fixture
def cwn():
    '''Returns a Coworking Network'''
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
    return CoworkingNetwork(cliques=collectors)

def test_cwn_nodes_count_attribute(cwn):
    '''Nodes count attribute keeps the number of times a collector appears in the dataset'''
    assert cwn.nodes['col1'].get('count')==8
    assert cwn.nodes['col4'].get('count')==4
    assert cwn.nodes['col5'].get('count')==2   
    assert cwn.nodes['col9'].get('count')==2
    
def test_cwn_isolated_nodes_included(cwn):
    '''Nodes without connections (isolated) are also included in the network'''
    assert 'col5' in cwn.nodes()
    
def test_cwn_count_attribute(cwn):
    pass
    
def test_initialize_cwn(cwn):
    '''Initialization'''
    assert cwn
    
def test_cwn_subgraphs_works(cwn):
    '''Derives subgraphs from connected components in the cwn model'''
    assert all( isinstance(sg,networkx.Graph) for sg in networkx.connected_component_subgraphs(cwn) )


# Execute tests above on script run
if __name__ == '__main__':
    pytest.main(['-v'])

#%%
g=cwn()
g.nodes['col1']
    
