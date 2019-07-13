import logging, re
from splice_list import SpliceList, Region

logger = logging.getLogger('trcls')

class Transcript:
    r"""
    A read transcript, e.g. a read pair.

    Attributes
    ----------
    segments : list
        ``list`` of ``Segment``\ s belonging to this transcript.

    Methods
    -------
    annotate(self, annotations, junction_tolerance)
        Annotates each ``Segment`` of this transcript according to the results of
        ``annotations.get_annotations``.
    """

    def __init__(self, sam_alignments, skip_tolerance, map_tolerance):
        r"""
        Parameters
        ----------
        sam_alignments : list
            ``list`` of SAM alignment lines representing the ``Segments`` that make
            up this transcript.
        skip_tolerance : int
            The minimum number of reference skips e.g. deletions, soft clips;
            for two regions at either end of the skip to be considered as
            separate exons. Alternatively, the maximum reference skip distance
            allowed in an exon.
        map_tolerance : int
            The minimum number of nucleotides mapped in a row before they are
            considered an exon.
        """
        self.segments = []

        for alignment in sam_alignments:
            try:
                self.segments.append(
                        Segment(alignment, skip_tolerance, map_tolerance))
            except NoMappingError:
                logger.warning(
                    '{} with flag {} has no mapping information. Skipping'
                    .format(*alignment.split('\t')[:2]))

        if self.segments == []:
            raise NoMappedSegmentsError

        # Set Junction complements. Junction complements are back-to-back exon
        # junctions in the transcript indicating that the transcript is a mature
        # mRNA spanning both exons on each side of that junction.
        c_lists = [segment.splice_list.components for segment in self.segments]
        components = [c for c_list in c_lists for c in c_list]
        for i in range(0, len(components)-1):
            if type(components[i]) == Region or type(components[i+1]) == Region:
                continue
            elif components[i].type != components[i+1].type and \
                    components[i].has_complement and \
                    components[i+1].has_complement:
                components[i].complement = components[i+1]
                components[i+1].complement = components[i]

        logger.debug('Transcript instantiated with\n\t<{}>'.format(
                '\n\t'.join(str(segment) for segment in self.segments)))

    def annotate(self, annotations, junction_tolerance):
        r"""
        Annotates each ``Segment`` of this transcript according to the results of
        ``annotations.get_annotations``.

        Parameters
        ----------
        annotations : Annotations
            Reference annotations.
        junction_tolerance : int
            Number of nucleotides which the transcript junctions are allowed to
            differ in position from a splice variant's junction such that the
            transcript can still be considered as being from the splice variant.
            This is particularly useful for transcripts aligned with non-splice
            aware aligners that have a tendency to slightly misalign transcripts
            near the splice junction.
        """
        represents = annotations.get_annotations(self, junction_tolerance)
        for segment in self.segments:
            segment.apply_annotation(represents)

    def __str__(self):
        return '\n'.join(segment.sam for segment in self.segments)

class Segment:
    r"""
    Represents a read fragment.

    Attributes
    ----------
    sam : str
        SAM alignment line representing this segment.
    splice_list : SpliceList
        Ordered list of regions and splice junctions for this segment.

    Methods
    -------
    apply_annotation(self, represents)
        Annotate this segment with the TR:Z optional field.
    """

    _CIGAR_MAPPING = ('M', 'D', 'N', '=', 'X')
    _CIGAR_MAPS_BOTH = ('M', '=', 'X')
    _CIGAR_MAPS_REF = ('D', 'N')

    def __init__(self, sam, skip_tolerance, map_tolerance):
        r"""
        Parameters
        ----------
        sam : str
            SAM alignment line representing this segment
        skip_tolerance : int
            The minimum number of reference skips e.g. deletions, soft clips;
            for two regions at either end of the skip to be considered as
            separate exons. Alternatively, the maximum reference skip distance
            allowed in an exon.
        map_tolerance : int
            The minimum number of nucleotides mapped in a row before they are
            considered an exon.
        """
        self.sam = sam
        region_list, set_junctions = self._get_region_list(
                sam, skip_tolerance, map_tolerance)
        if region_list == []:
            # the mapped regions exist but are all within map_tolerance,
            # so nonw are counted as a 'region'
            raise NoMappingError
        self.splice_list = SpliceList(
                self.sam.split('\t')[0], region_list, **set_junctions)

    def __str__(self):
        return '{} <{}>'.format(_get_qname(self.sam), str(self.splice_list))

    def apply_annotation(self, represents):
        r"""
        Annotate this segment with the TR:Z optional field. If there are no
        representations (i.e. an empty list is passed), the value of the TR:Z
        field is set as ``TR:Z:*``.

        e.g. ``TR:Z:NM_001456`` tags a segment belonging to transcript from
        NM_001456.
        e.g. ``TR:Z:NM_001456,NM_001110556`` tags a segment belonging to
        transcript from either NM_001456 or NM_001110556.
        e.g. ``TR:Z:*`` if no information available (failed to match any variant
        in annotation, include pre-mRNA)

        Parameters
        ----------
        represents : list
            ``list`` of ``str``, each string is an accession number of a splice
            variant OR 'pre-mRNA' if this transcript originates from pre-mRNA.
        """

        TR_tag = 'TR:Z:{}'.format(
                '*' if represents == [] else ','.join(represents))
        self.sam = '\t'.join([self.sam, TR_tag])

    @classmethod
    def _get_region_list(cls, sam_alignment, skip_tolerance, map_tolerance):
        r"""
        Get a list of regions which this segment spans.
        """
        position = int(sam_alignment.split('\t')[3])
        cigar_v = cls._expand_cigar(sam_alignment.split('\t')[5])

        # Remove ops mapping the transcript only (e.g. insertions into ref); but
        # retain soft_clips to determine if SpliceList should mark the left or
        # right or both ends as Junctions.
        cigar_v = filter(lambda op: op in cls._CIGAR_MAPPING + ('S',), cigar_v)
        cigar_v = ''.join(cigar_v)
        left_soft_clip = re.match(r'\AS+', cigar_v)
        should_set_left_junction = False if left_soft_clip == None \
                else len(left_soft_clip.group()) > skip_tolerance
        right_soft_clip = re.search(r'S+\Z', cigar_v)
        should_set_right_junction = False if right_soft_clip == None \
                else len(right_soft_clip.group()) > skip_tolerance
        # Once soft clips are checked, remove them
        cigar_v = filter(lambda op: not op == 'S', cigar_v)
        # Encode non-skips as True and skips as False
        cigar_v = map(
                lambda op: True if op in cls._CIGAR_MAPS_BOTH else False, cigar_v)
        # Re-represent as list<tuple<boolean, occurences>>
        cigar_buffer = []  # elem: tuple<value, occurences>
        times_skipped = 0
        for op in cigar_v:
            # Init step
            if cigar_buffer == []:
                cigar_buffer.append([op, 1])
            # If this op is the same as the last
            elif op == cigar_buffer[-1][0]:
                cigar_buffer[-1][1] += 1
            # This op is not the same as the last
            else:
                cigar_buffer.append([op, 1])
        # Transform skips to non-skips according to skip_tolerance
        cigar_buffer = map(
                lambda pair: (True, pair[1]) if \
                        pair[0] == False and pair[1] <= skip_tolerance else \
                        pair,
                cigar_buffer)
        # Merge adjacent non-skips
        cigar_buffer2 = []
        for pair in cigar_buffer:
            # Init step
            if cigar_buffer2 == []:
                cigar_buffer2.append(pair)
            # Should merge
            elif pair[0] == cigar_buffer2[-1][0]:
                cigar_buffer2[-1][1] += pair[1]
            # Shouldn't merge
            else:
                cigar_buffer2.append(pair)
        # Translate skip and non-skip as regions. Follow the GTF indexing:
        # 1-based, inclusive on both sides.
        regions = []
        for region in cigar_buffer2:
            if region[0] == True:
                regions.append([position, position+region[1]-1])
                position += region[1]
            else:
                position += region[1]

        # remove exons which span < map_tolerance
        within_tolerance = []
        for region in regions:
            if not region[1] - region[0] + 1 < map_tolerance:
                within_tolerance.append(region)

        return (within_tolerance, {
                        'set_left_junction': should_set_left_junction,
                        'set_right_junction': should_set_right_junction})

    @staticmethod
    def _expand_cigar(cigar):
        r"""
        Expands a CIGAR string such that the length of the string is the number
        of CIGAR operations, and each character at its index is equal to the
        operation at that CIGAR 'step'.
        """
        if cigar == '*':
            raise NoMappingError
        expanded = ''
        number_of_operations = None  # also a flag for the parse state
        for char in cigar:
            # the start of a number
            if char.isnumeric() and number_of_operations == None:
                number_of_operations = char
            # the continuation of a number
            elif char.isnumeric():
                number_of_operations += char
            # the naming of an operation
            else:
                expanded += int(number_of_operations) * char
                number_of_operations = None
        return expanded


class NoMappingError(Exception):
    r"""
    Exception to be raised when an alignment does not contain any mapping
    information. That is, the CIGAR field == '*'.
    """
    pass


class NoMappedSegmentsError(NoMappingError):
    r"""
    Exception to be raised when an alignment group (i.e. transcript) does not
    contain any mapping information. E.g. a paired read has both pairs unmapped.
    """
    pass


def get_transcripts(alignments, skip_tolerance, map_tolerance):
    r"""
    Get the ``Transcript``\ s for many alignments.

    Parameters
    ----------
    alignments : sequence
        Sequence iterable of SAM alignment lines. Each line should be stripped
        (i.e. no newline character at the end).
    skip_tolerance : int
        The minimum number of reference skips e.g. deletions, soft clips; for
        two regions at either end of the skip to be considered as separate
        exons. Alternatively, the maximum reference skip distance allowed in an
        exon.
    map_tolerance : int
        The minimum number of nucleotides mapped in a row before they are
        considered an exon.

    Returns
    -------
    map
        Iterable map object of ``Transcript``\ s.
    """
    transcript_groups = _get_read_groups(alignments)
    transcripts = []
    for group in transcript_groups:
        try:
            transcripts.append(Transcript(group, skip_tolerance, map_tolerance))
        except NoMappedSegmentsError:
            pass
    return transcripts

def _get_read_groups(alignments):
    r"""
    Group alignments according to read groups (most commonly read pairs,
    or a single-member group for unpaired reads).

    Parameters
    ----------
    alignments : sequence
        Sequence iterable of SAM alignment lines.

    Returns
    -------
    list
        ``list`` of ``tuple``, where each tuple contains a sequence of strings
        belonging to the same read group; each string is a full SAM alignment
        line.
    """

    alignments = list(map(str.strip, alignments))

    is_grouped = lambda line: _get_flags(line) & 0x1 == 0x1
    is_ungrouped = lambda line: not is_grouped(line)
    is_first_segment = lambda line: _get_flags(line) & 0x40 == 0x40
    is_last_segment = lambda line: _get_flags(line) & 0x80 == 0x80

    # Extract tthe alignments which are not in read groups
    alignments_ungrouped = list(filter(is_ungrouped, alignments))

    # Extract he alignments which are in read groups
    alignments_grouped = list(filter(is_grouped, alignments))
    grouped_first = list(filter(is_first_segment, alignments_grouped))
    grouped_last = list(filter(is_last_segment, alignments_grouped))
    grouped_middle = list(filter(
            lambda seg: not is_first_segment(seg) and  not is_last_segment(seg),
            alignments_grouped))
    num_grouped = {'total_segments': len(alignments_grouped), 'has_group': 0}

    # Get read groups by starting from the first segment in the group, following
    # the RNEXT tag, until a last segment in the group in found.
    grouped_orphans = []
    read_groups = []
    read_group = []
    for alignment in grouped_first:
        read_group = [alignment]
        rnext = _get_rnext(alignment)
        rnext = _get_qname(alignment) if rnext == '=' \
                else None if rnext == '*' else rnext

        # Consider grouped_middles
        i = 0
        while i < len(grouped_middle) and rnext != None:
            consider = grouped_middle[i]
            consider_qname = _get_qname(consider)
            # The i'th read is the next segment.
            if consider_qname == rnext:
                read_group.append(consider)
                rnext = _get_rnext(consider)
                rnext = _get_qname(consider) if rnext == '=' \
                        else None if rnext == '*' else rnext
                grouped_middle.pop(i)
            # The i'th read is not the next segment.
            else:
                i += 1

        # Consider grouped_last
        for i, consider in enumerate(grouped_last):
            consider_qname = _get_qname(consider)
            if consider_qname == rnext:
                read_group.append(consider)
                grouped_last.pop(i)
                break

        # If no next segment can be found from a first segment, consider it
        # orphaned
        if len(read_group) == 1:
            grouped_orphans.append(read_group[0])
            continue
        elif not is_last_segment(read_group[-1]):
            logger.warning(
                    'Read group with first segment {} does not end with a '
                    'last segment.'.format(_get_qname(read_group[0])))

        read_groups.append(read_group)
        num_grouped['has_group'] += len(read_group)

    # Reads that were marked as group, but not belonging to any read group.
    # TODO: find groups within the leftovers.
    grouped_orphans += grouped_middle + grouped_last

    logger.info(
            'Found {} total reads. {} were ungrouped. {} were grouped into '
            '{} read groups. {} were orphaned reads.'.format(
                    len(alignments), len(alignments_ungrouped),
                    num_grouped['has_group'], len(read_groups),
                    len(grouped_orphans)))

    return read_groups + [[alignment] for alignment in alignments_ungrouped] + \
            [[alignment] for alignment in grouped_orphans]

def _get_flags(alignment):
    r"""
    Given a SAM alignment line, get the bitwise FLAGs as int
    """
    return int(alignment.split('\t')[1])

def _get_qname(alignment):
    r"""
    Given a SAM alignment line, get the QNAME as str. Checks if alignment ==
    None to facillitate the iteration over other_segments in
    get_alignment_groups.
    """
    if alignment == None:
        return None
    return alignment.split('\t')[0]

def _get_rnext(alignment):
    r"""
    Given a SAM alignment line, get the RNEXT field as str
    """
    return alignment.split('\t')[6]
