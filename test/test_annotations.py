import os, sys
import pytest

sys.path.append(os.path.join(sys.path[0], '../src'))

from test_annotations_cases import *
from splice_list import SpliceList
from annotations import Annotations
from transcript import Transcript

regions_variants = [[NM_001110556_exons, NM_001456_exons, pre_mRNA_whole]]

@pytest.fixture
def flna_annotations():
    with open('test/FLNA.gtf') as gtf_file:
        return Annotations(gtf_file)

@pytest.mark.parametrize('regions_variants', regions_variants)
def test_annotations(flna_annotations, regions_variants):
    splice_lists_variants = [
            SpliceList('', regions, set_left_junction=True,
                set_right_junction=True)
            for regions in regions_variants]

    for annotation, variant in zip(
            flna_annotations.variants, splice_lists_variants):
        assert annotation.contains(variant, 0)
        assert variant.contains(annotation, 0)

@pytest.mark.parametrize('args', pre_mRNA_only)
def test_pre_mRNA_only(flna_annotations, args):
    transcript = Transcript(*args['transcript_args'])
    assert flna_annotations.get_annotations(
            transcript, args['junction_tolerance']) == ['pre-mRNA']

@pytest.mark.parametrize('args', mature_mRNA_only)
def test_mature_mRNA_only(flna_annotations, args):
    transcript = Transcript(*args['transcript_args'])
    assert flna_annotations.get_annotations(
            transcript, args['junction_tolerance']) == ['NM_001110556', 'NM_001456']

@pytest.mark.parametrize('args', NM_001110556_only)
def test_NM_001110556_only(flna_annotations, args):
    transcript = Transcript(*args['transcript_args'])
    assert flna_annotations.get_annotations(
            transcript, args['junction_tolerance']) == ['NM_001110556']

@pytest.mark.parametrize('args', NM_001456_only)
def test_NM_001456_only(flna_annotations, args):
    transcript = Transcript(*args['transcript_args'])
    assert flna_annotations.get_annotations(
            transcript, args['junction_tolerance']) == ['NM_001456']

def test_completely_ambiguous(flna_annotations):
    pass

def test_no_matches(flna_annotations):
    pass
