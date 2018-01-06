#!/usr/bin/env python

"""
Names Cleaning module
"""
from . import NamesMap
import json, re, string, unicodedata
from collections import Counter
from warnings import warn





# =================
# Names Atomization
# -----------------

def namesFromString( namesStr, delim=';', unique=False, preserveOrder=False ):
    """
    Atomizes names in a names string using specified delimiters. 
    
    Parameters
    ----------
    namesStr : str
        Names string to be atomized.
        
    delim : str or list, default ';'
        Delimiter that is used to separate names in the names string. 
        If a list of delimiters is passed they will be concatenated
        into a regular expression for splitting.
    
    unique : bool, default False
        If set to true, returned values are unique. 
    
    preserveOrder : bool, default False
        If True the order in which names appear in the string is 
        guaranteed to be preserved in the list (slower).
        
    Returns
    -------
        A list with the names extracted from the string.
    """
    if type(delim)==str:
        namesSplit = namesStr.split(delim)
        
    elif type(delim)==list:
        delim = [ '\|' if i=='|' else i for i in delim ] 
        namesSplit = re.split( '|'.join( c for c in delim ), namesStr)
    
    namesList = [ n for n in [ name.strip() for name in namesSplit ] if n!='' ]
    
    if unique:
        if not preserveOrder: return list(set(namesList))
        else:
            namesCounts = dict( (n,0) for n in namesList )
            unique_namesList = []
            for n in namesList:
                if namesCounts[n]: continue
                namesCounts[n]+=1
                unique_namesList.append(n)
        return unique_namesList   
    
    return namesList


def atomizeNames( col, operation=None, replaces=None ):
    """
    Applies an atomization operation on a names column, which must be a pandas Series. 
    The atomized names at each row are stored as a list.
    
    Parameters
    ----------
    col : pandas.Series
        Names column to be atomized.
        
    operation : function
        The atomization operation to be applied to the names column
        
    replaces:
        A list of 2-tuples (srclst, tgt), where srclst is a list of names to be replaced by tgt.
        The element tgt can be either a string or a function which results in a string.
        
    Returns
    -------
        A pandas Series with lists of atomized names.
    """
    if replaces is not None:
        replacesDict = dict( (src, tgt(src)) if callable(tgt) else (src,tgt) for (srclst, tgt) in replaces for src in srclst )
        col = col.replace( replacesDict )
        
    col_atomized = col.apply( operation )
    return col_atomized


def getNamesList( col, with_counts=False, orderBy=None ):
    """
    Gets a list of names from an atomized names column.
    
    Parameters
    ----------
    col : pandas.Series
        Atomized names column from which to retrieve names.
    
    with_counts : bool
        If set to True the result includes the number of records by each collector.
    
    orderBy : str
        If some rule is specified, the resulting list is sorted. Rules can be either
        to sort alphabetically ('alphabetic') or by the number of records by each
        collector ('counts').   
    """
    
    if orderBy not in [ None, "alphabetic", "counts"]:
        raise ValueError("Invalid argument for 'orderBy': {}".format(orderBy))
    
    if with_counts or orderBy=="counts":
        l = [ (n,c) for (n,c) in Counter( n for nlst in col for n in nlst ).items() ] 
        if orderBy=="alphabetic":
            return sorted( l, key=lambda x: x[0] )
        elif orderBy=="counts":
            sorted_l = sorted( l, key=lambda x: x[1], reverse=True )
            if with_counts:
                return sorted_l
            else:
                return [ n for (n,c) in sorted_l ]
        else:
            return l
                
    else:
        if orderBy=="alphabetic":
            return sorted(list(set( n for nlst in col for n in nlst )))
        else:
            return list(set( n for nlst in col for n in nlst ))


# ===================
# Names normalization
# -------------------

def normalize(name, normalizationForm='NFKD'):
    """
    A simple normalization function. A name is split on commas, periods and spaces are removed and it is then set to a unicode normalization form.
    
    Parameters
    ----------
    name : str
        Name to be normalized.
    
    normalizationForm: str, default 'NFKD'
        Unicode normalization form to apply during the process.
        
    Examples
    --------
    "João da Silva" -> "joaodasilva"
    "Leite, A.M." -> "leite,am"
    "J. R. Souza" -> "jrsouza"
    """
    name = name.lower() # to lowecase
    name = name.replace('.','') # remove periods
    name_ls = tuple( part.strip() for part in name.split(',') ) # split and strip names into tuples

    normalize = lambda s: ''.join( x for x in unicodedata.normalize(normalizationForm, s) if x in string.ascii_letters ) # remove accents
    name_ls = tuple( normalize(name) for name in name_ls )
    
    return ','.join(name_ls)


# ------------------
# Names Map
# ------------------


            

def read_NamesMap_fromJson(filepath, normalizationFunc=None):
    """
    Creates a NamesMap instance from a json file containing both a primitive-to-normalized names map
    and a remapping index. The json object must have both attributes `map_prim_norm` and 
    `_remappingIndex`, which stores data used to instance NamesMap class.
    
    Parameters
    ----------
    filepath : str
        Path to the json file containing the map
        
    normalizationFunc : function
        A normalization function to be passed to the NamesMap constructor. If it is 
        not set a warning is issued, as the NamesMap will not be assigned to any
        normalization rule.
    """
    if normalizationFunc is None:
        warn("A names map was created without a normalization function!")
        
    with open(filepath,'r') as f:
        d = json.load(f)
        nm = NamesMap( names=None, normalizationFunc=normalizationFunc, 
                       _map_prim_norm=d['_map_prim_norm'], 
                       _remappingIndex=d['_remappingIndex'])
    
    return nm



# ==============
# Names indexing
# --------------

def getNamesIndexes( df, atomizedNamesCol, namesMap=None ):
    # split_names_col is a column with names already split
    namesIndexes = dict( (name,[]) for name in namesMap.values() )
    for i,names in df[atomizedNamesCol].iteritems():
        for name in names:
            if namesMap is not None:
                try:
                    namesIndexes[namesMap[name]].append(i)
                except KeyError:
                    pass
            else:
                namesIndexes[name].append(i)
            
    return namesIndexes
