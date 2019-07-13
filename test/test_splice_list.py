import os, sys
import pytest

sys.path.append(os.path.join(sys.path[0], '../src'))

from test_splice_list_cases import *
from splice_list import SpliceList, Region, Junction, OverlappingExonsError

@pytest.mark.parametrize('identifier,exons', exons_list)
def test_splice_list(identifier, exons):
    splice_list = SpliceList(identifier, exons)
    assert splice_list.identifier == identifier
    assert len(splice_list.components) == 3 * len(exons) - 2
    assert len(splice_list.regions) == len(exons)
    assert type(splice_list.components[0]) == Region
    assert type(splice_list.components[-1]) == Region

@pytest.mark.parametrize('exon', exon_list)
def test_region(exon):
    region = Region(*exon)
    assert region.start == exon[0]
    assert region.stop == exon[1]
    assert str(region) == 'region {}-{}'.format(*exon)

@pytest.mark.parametrize('identifier,exons', overlapping_exons_list)
def test_splice_list_raises_overlapping_exons_error(identifier, exons):
    with pytest.raises(OverlappingExonsError):
        splice_list = SpliceList(identifier, exons)

@pytest.mark.parametrize('position,_type', junction_list)
def test_junction(position, _type):
    junction = Junction(position, _type)
    assert junction.position == position
    assert junction.type == _type
    assert str(junction) == 'junction ({}) at {}'.format(
            'start' if _type == Junction.TYPE_START else 'end', position)

@pytest.mark.parametrize('exons1,exons2,match_tolerance', matching_exons)
def test_splice_list_contains_True(exons1, exons2, match_tolerance):
    splice_list1 = SpliceList('test', exons1)
    splice_list2 = SpliceList('test', exons2)
    assert splice_list1.contains(splice_list2, match_tolerance) == True

@pytest.mark.parametrize('exons1,exons2,match_tolerance', non_matching_exons)
def test_splice_list_contains_True(exons1, exons2, match_tolerance):
    splice_list1 = SpliceList('test', exons1)
    splice_list2 = SpliceList('test', exons2)
    assert splice_list1.contains(splice_list2, match_tolerance) == False
