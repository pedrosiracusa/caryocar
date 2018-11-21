Caryocar
========

Caryocar is a Python package for building Species-Collector Networks (SCNs) and Collector CoWorking Networks(CWNs) models from species occurrence data, as introduced in my [MSc thesis](https://lncc-netsci.github.io/pedrocs/).
SCNs and CWNs extend the social network analytics and can be used for understanding the social structure behind biological collections.
This package is built on top of [NetworkX](https://github.com/networkx/networkx).


### Supporting documents
* New perspectives on analyzing data from biological collections based on social network analytics [[MSc thesis](https://tede.lncc.br/handle/tede/279)].
* On the social structure behind biological collections [[Preprint](https://www.biorxiv.org/content/early/2018/06/08/341297)].
* Package documentation coming soon...


Example Usage
-------------
Create a Species-Collector Network (SCN) from a list of collectors and species:

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



Create a Collector CoWorking Network (CWN) from a list of collector cliques:

		>>> collectors = [ ['a','b','c'], ['d','e'], ['a','c'] ]
    >>> cwn = CoworkingNetwork(cliques=collectors)
    
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


Install
-------

This package is still experimental, and should ideally be run from a conda virtual environment, which is specified in the `environment.yml` file. In order to create the virtual environment clone this repository, make sure you have [conda](https://conda.io/docs/) installed and use one of the following commands, from the root of the repository:

    $ conda env create -f environment.yml

Then you should activate it with the following command:

* On Linux:


    $ source activate caryocar


* On Windows:


    $ activate caryocar

