import pytest
from caryocar.cleaning import NamesMap

# ===================
# test NamesMap class
# ===================
@pytest.fixture
def nm():
    names = ["name"+str(i) for i in range(1,10)]
    upper = lambda x: x.upper()
    return NamesMap(names=names, normalizationFunc=upper)

@pytest.fixture
def nm_remapped():
    names = ["name"+str(i) for i in range(1,10)]
    remap = { 'NAME1':'NAME2',
              'NAME2':'NAME3',
              'NAME3':'name_3',
              'NAME4':'NAME_4',
              'NAME5':'name5',
              'NAME6':'NAME5' }
    upper = lambda x: x.upper()
    
    return NamesMap(names=names, normalizationFunc=upper, remappingIndex=remap)


def test_namesmap_get_normalized_names(nm):
    '''Tests the retrieval of normalized names map'''
    assert all( nm._normalizationFunc(orig)==norm for (orig,norm) in nm.getMap().items() )

@pytest.mark.parametrize("origName,remappedName",[
        ('name1','name_3'),
        ('name4','NAME_4')
        ])
def test_namesmap_remapped_getMap_remap(nm_remapped,origName,remappedName):
    '''The names map getMap method performs remapping if the argument is true'''
    assert nm_remapped.getMap(remap=True)[origName]==remappedName
    assert nm_remapped.getMap(remap=False)[origName]==nm_remapped._normalizationFunc(origName)
    
@pytest.mark.parametrize("remap",[
        { 'NAME1':'NAME2',
          'NAME2':'NAME1' },
        { 'NAME1':'NAME2',
          'NAME2':'NAME3',
          'NAME3':'NAME1' },
        { 'NAME1':'NAME2',
          'NAME2':'NAME3',
          'NAME3':'NAME4',
          'NAME4':'NAME2' } ])
def test_namesmap_remapped_getMap_loopbacks_raise_runtimeerror(remap):
    '''The names map getMap method raises a runtime error if a loopback in the remapping index is detected and remap argument is set True'''
    names=["name"+str(i) for i in range(1,10)]
    upper = lambda x: x.upper()
    nm=NamesMap(names=names,normalizationFunc=upper,remappingIndex=remap)
    with pytest.raises(RuntimeError):
        nm.getMap(remap=True)
    assert nm.getMap(remap=False)

@pytest.mark.parametrize("remap,endpoint",[
        ({ 'NAME1':'NAME2',
          'NAME2':'NAME1' },'NAME2'),
        ({ 'NAME1':'NAME2',
          'NAME2':'NAME3',
          'NAME3':'NAME1' },'NAME2'),
        ({ 'NAME1':'NAME2',
          'NAME2':'NAME3',
          'NAME3':'NAME4',
          'NAME4':'NAME2' },'NAME2') ])    
def test_namesmap_remapped_loopbacks_setting_endpoints(remap,endpoint):
    '''Setting endpoints fixes loopbacks in names maps'''
    names=["name"+str(i) for i in range(1,10)]
    upper=lambda x: x.upper()
    nm=NamesMap(names=names,normalizationFunc=upper,remappingIndex=remap)
    nm.setEndpoint(endpoint)
    assert nm.getMap()
    assert any( upper(u)==v for u,v in nm.getMap().items() if v==endpoint ) # the endpoint should have its non-normalized form as a source
    
# Execute tests above on script run
if __name__ == '__main__':
    pytest.main(['-v', __file__ ])