#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Names Atomizer module
"""

import json
from collections import Counter

class NamesAtomizer:
    """
    The NamesAtomizer is built with an atomizing operation to be defined as the instance's default and an optional list with names to be replaced. Names to be replaced must be passed in a list of tuples, in any of the following ways:
            
    >> rep = [('n1', 'correct_n1'), ('n2', 'correct_n2')]
    or
    >> rep = [(['n1','n1_2], 'correct_n1'), (['n2'], 'correct_n2')]
    or
    >> expr1 = lambda x: x.replace(';', '_')
    >> expr2 = lambda x: x.replace('&', '_')
    >> rep = [(['n1;1', 'n1;2'], expr1), (['n2&1', 'n2&2'], expr2)]
    
    Note that if you pass an expression as the second item of the tuple this expression must evaluate in a string!
    
    Class methods
    -------------
    .atomize
    .addReplaces
    .writeReplaces
    .readReplaces
    .getCachedNames
    """
        
    def __init__(self, atomizeOp, replaces=None):
        """
        Initialization of NamesAtomizer class.
        Parameters
        ----------
        atomizeOp : function
        replaces: list of tuples
        """
        self._replaces = self._buildReplaces(replaces)
        self._operation = atomizeOp
        self._cache = None
        
    def _buildReplaces(self, replacesList):
        """
        Builds a replaces dict from a list. The input list must contain tuples in which the first element is a list of names strings that must be replaced by the string in the tuple's second element. The second element can alternatively be an expression that results in a names string.
        """
        res = {}
        if replacesList is None:
            return res
        
        for srcs,tgt in replacesList:
            if isinstance(srcs, str):
                src = srcs
                res.update( {src:tgt(src)} if callable(tgt) else {src:tgt})
                
            elif isinstance(srcs, (list,tuple,set)):
                for src in srcs:
                    res.update( {src:tgt(src)} if callable(tgt) else {src:tgt} )

            else:
                raise ValueError("Invalid value '{0}' in '({0},{1})'. Must be either string or iterable".format(srcs, tgt))
                
        return res                 
    
    def atomize(self, col, operation=None, withReplacing=True, cacheResult=True):
        """
        This method takes a column with names strings and atomizes them
        
        Parameters
        ----------
        
        col : pandas.Series
            A column containing names strings to be atomized.   
        
        operation : function
            If an operation is passed in it is used to atomize the column instead
            of the instance's default operation.
        
        withReplacing : bool, default True
            If set to True some names replacing is performed before atomization. 
        
        cacheResult : bool, default True
            If set to True the resulting series is cached for later use.
        """
        if operation is None:
            operation = self._operation
        
        if withReplacing:
            col=col.replace(self._replaces)
            
        atomizedCol = col.apply(operation)
        if cacheResult:
            self._cache = (col, atomizedCol)
        return atomizedCol
    
    def addReplaces(self, replacesList):
        replacesDict = self._buildReplaces(replacesList)
        self._replaces.update(replacesDict)
    
    def write_replaces(self, filename):
        """
        Writes replaces to a json file
        """
        with open(filename,'w') as f:
            d = {'_replaces':self._replaces}
            json.dump(d, f, sort_keys=True, indent=4, ensure_ascii=False)
            
    def read_replaces(self, filepath, update=True):
        """
        Reads replaces from a json file
        """
        with open(filepath, 'r') as f:
            data = json.load(f)
            if update:
                self._replaces.update( data['_replaces'] )
            else:
                self._replaces = data['_replaces']
        
    def getCachedNames(self, namesToFilter=['et al.'], sortingExp=lambda x: [len(x[0]),-x[2]]):
        """
        This method uses data in the instance's cache. Returns atomized names from the instance's cache. Names are associated to their original namestring as well as the number of records they appear in the dataset. The result is structured as a 3-tuple, with elements in the same order stated above.
        
        Parameters
        ----------
        
        namesToFilter : list
            Names that should be ignored by the method. By default it ignores 'et al.'.
        
        sortingExp : function
            An expression to be passed as key to sort the final result.
            
        Returns
        -------
        
        A 3-tuple (u,v,w) where:
          u = atomized name;
          v = original name string that was used to atomize names;
          w = count of the total occurrences of an atomized name in the dataset.
        """
        c = self._cache
        l = [ (n,nstr) for nstr,norm in zip(c[0],c[1]) for n in norm if n not in namesToFilter ]
        ctr = Counter(i[0] for i in l)
        return sorted([ (u,v,ctr[u]) for u,v in set(l) ],key=sortingExp)