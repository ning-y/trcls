import os, sys
import pytest

sys.path.append(os.path.join(sys.path[0], '../trcls'))

from test_transcript_cases import *
from transcript import Transcript, Segment, _get_read_groups, \
        NoMappingError, NoMappedSegmentsError
from cli import DEFAULT_JUNCTION, DEFAULT_MAP, DEFAULT_SKIP

@pytest.mark.parametrize('cigar,cigar_expanded', cigar_and_cigar_expanded)
def test_Segment_expand_cigar(cigar, cigar_expanded):
    assert Segment._expand_cigar(cigar) == cigar_expanded

@pytest.mark.parametrize('sam,regions,skip_tolerance,map_tolerance', sam_and_regions)
def test_Segment_get_region_list(sam, regions, skip_tolerance, map_tolerance):
    assert Segment._get_region_list(sam, skip_tolerance, map_tolerance) == regions

@pytest.mark.parametrize('alignments,groups', transcript_groups)
def test_get_read_groups(alignments, groups):
    assert _get_read_groups(alignments) == groups

@pytest.mark.parametrize('alignment', no_mapping)
def test_raises_no_mapping_error(alignment):
    with pytest.raises(NoMappingError):
        Segment(alignment, DEFAULT_SKIP, DEFAULT_MAP)

@pytest.mark.parametrize('alignments', no_mapped_segments)
def test_raises_no_mapped_segments(alignments):
    with pytest.raises(NoMappedSegmentsError):
        Transcript(alignments, DEFAULT_SKIP, DEFAULT_MAP)
