#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Names Map module
"""

import json
from copy import deepcopy
from warnings import warn
from collections import Counter


class NamesMap:
    """
    The class which describes NamesMap objects. Names maps store both name primitives and normalized names. Primitives are the "original names", as they're given as input to the class constructor. When the class is instanced each name primitive is mapped to its normalized form through a normalization function. Normalized maps can then be remapped to other names by following a remapping index.
    
    Class methods
    -------------
    .getInconsistencies
    .getMap
    .addNames
    .remap
    .setEndpoint
    .write_toJson
    """
    
    def __init__(self, names, normalizationFunc, remappingIndex=None, *args, **kwargs):
        """
        The NamesMap constructor
        
        Parameters
        ----------
        names : list
            A list with names to be normalized
        
        normalizationFunc : function
            A function or expression for names normalization
            
        remappingIndex : dict
            A dictionary with mapped names to initialize the instance's remapping index
        """
        self._normalizationFunc = normalizationFunc
        
        load_map_prim_norm = kwargs.get('_map_prim_norm',None)
        load_remappingIndex = kwargs.get('_remappingIndex',None)
        
        normNamesDict = lambda nlst: dict( (n,self._normalizationFunc(n)) for n in nlst )
        self._map_prim_norm = normNamesDict(names) if load_map_prim_norm is None else load_map_prim_norm 
        self._remappingIndex = remappingIndex if load_remappingIndex is None else load_remappingIndex
    

    def _getRef(self,n):
        """
        Follows all chained references for a name in the remapping index
        
        Parameters
        ----------
        n : str
            The name to be de-referenced
        """
        start=n
        chain=[]
        remappingKeys = self._remappingIndex.keys()
        while n in remappingKeys:
            chain.append(n)
            n = self._remappingIndex[n]
            if n in chain:
                chain.append(n)
                raise RuntimeError("Loopback detected", start, chain)
        return n
    
    
    def _get_loopback_inconsistencies(self):
        """
        Detects loopbacks in in mapping chains
        """
        inconsistencies = {}
        for k in self._remappingIndex.keys():
            try:
                self._getRef(k)
                
            except RuntimeError as e:
                inconsistencies['mes'] = inconsistencies.get('mes',[]) + [e.args[0]]
                inconsistencies['key'] = inconsistencies.get('key',[]) + [e.args[1]]
                inconsistencies['chain'] = inconsistencies.get('chain',[]) + [e.args[2]]
        
        if len(inconsistencies)==0:
            return None
        else:
            return inconsistencies
    
    
    def _remove_selfloops(self):
        keys_to_remove = [ k for k in self._remappingIndex.keys() if k==self._remappingIndex[k]]       
        for k in keys_to_remove:
            self._remappingIndex.pop(k)
                  
            
    def getInconsistencies(self, prettyPrint=True):
        d = {}
        d['loopback_inconsistencies'] = self._get_loopback_inconsistencies()
        
        if any( True if v is not None else False for v in d.values()  ):
            if prettyPrint:
                mes = "INCONSISTENCIES\n===============\n"

                # loopback inconsistencies
                if d['loopback_inconsistencies'] is not None:
                    mes += "Loopback Inconsistencies\n"
                    data = list(zip( *d['loopback_inconsistencies'].values() ))
                    dataStr = lambda t: "  > {}: Starting from key '{}' got chain {}\n".format(*t)
                    mes += ''.join( dataStr(t) for t in data )
                    mes += '---------------'

                return mes
            
            return d
        
        return None
                
        
    def getMap(self, remap=True):
        """
        Returns a COPY of the names map.
        
        Parameters
        ----------
        remap : bool, default True
            If set to True, the names map is buit by first de-referencing
            remaps in the remapping index. Otherwise all remaps will not
            be considered for building the names map.
        """
        res = deepcopy(self._map_prim_norm)
        if remap and self._remappingIndex is not None:
            for s,t in self._map_prim_norm.items():
                try:
                    res[s] = self._getRef(t)
                except RuntimeError as e:
                    raise(e)
        return res
    
    def addNames(self, names, normalizationFunc=None, updateExistingKeys=False):
        """
        Updates the names map using a list of primitive names, which are stored as references
        to their normalized forms.
        
        Parameters
        ----------
        names : list
            List with names to be inserted or updated in the map. They are stored as primitives,
            mapping to their normalized forms.
            
        normalizationFunc : function
            By default the instance's normalization function is used to normalize the new names 
            primitives. If an alternative expression is passed in it will be used to perform the
            normalization instead.
            
        updateExistingKeys : bool, default False
            By default only names that still do not exist as keys in the primitives-normalized 
            names map are normalized and updated. If set to True all input names will be 
            updated in the names map.  
        """
        normFunc = self._normalizationFunc if normalizationFunc is None else normalizationFunc
        d = dict( ( n,normFunc(n) ) for n in names )
        
        if not updateExistingKeys:
            d = dict( (k,v) for k,v in d.items() if k not in self._map_prim_norm.keys() )
        
        self._map_prim_norm.update(d) 
            
    
    def remap(self, remaps, fromScratch=False):
        """
        Updates the remapping dictionary using a list of tuples as input.
        
        Parameters
        ----------
        remaps : list of tuples
            Remaps values from tuples (s,t), where a normalized name s remaps to a
            normalized name t.
        
        fromScratch : bool
            If set to True the remapping dict becomes the one passed in. All other previous
            remaps are discarded.
            
        Note
        ----
        If the list of tuples passed in contains duplicated keys a warning is issued, and the
        latest (key,value) pair is the one which will persist.
        """        
        # check for duplicated keys
        duplicatedKeys = [ s for s,cnts in Counter( s for s,t in remaps ).items() if cnts>1 ]
        if len(duplicatedKeys)>0: 
            warningMsg = "Some keys from input are duplicated: {}.".format(str(duplicatedKeys))
            warn(warningMsg)
        
        # update remapping index
        if fromScratch: self._remappingIndex=None
        if self._remappingIndex is None: self._remappingIndex={}
        
        for s,t in remaps:
            self._remappingIndex[s] = t
        
        self._remove_selfloops()
        return self.getInconsistencies()
    
    
    def setEndpoint(self, key):
        """
        This method sets a name as the endpoint of a chain. This method is used for resolving loopbacks
        in the remapping chain. The name is set to be the latest reference in the chain, and therefore
        does not map to any other name.
        
        Parameter
        ---------
        key : str
            The name to be set as the latest reference.
        """
        return self._remappingIndex.pop(key)
    
    
    def write_toJson(self, filename, flatten=False):
        """
        Creates a json file to store a NamesMap's primitive-to-normalized names map and remapping index.
        Data is stored as json object arguments `_map_prim_norm` and `_remappingIndex`.
        
        Parameters
        ----------
        filename : str
            Path to the file to be created.
        
        flatten : bool, default False
            If set to true, all remappings are consolidated into the `_map_prim_norm` map. In other words,
            the remapping index is used to assign each name primitive to its final reference. All references
            are then removed from the remapping index.
        
        """
        json_dict = dict([ ('_map_prim_norm', self.getMap() if flatten else self._map_prim_norm),
                           ('_remappingIndex', {} if flatten else self._remappingIndex) ])
        
        with open(filename, 'w') as output_file:
            json.dump( json_dict, output_file, sort_keys=True, indent=4, ensure_ascii=False)           