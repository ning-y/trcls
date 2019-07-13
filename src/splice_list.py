import logging

logger = logging.getLogger('trcls')

class SpliceList:
    r"""
    Ordered list of expressed regions and splice junctions representing a splice
    variant or precursor mRNA.

    Attributes
    ----------
    identifier : str
        Accession number if this ``SpliceList`` represents a splice variant,
        otherwise 'pre-mRNA' for precursor mRNAs.
    components : list
        Ordered ``list`` of ``Region``\ s and ``Junction``\ s.
    regions : list
        Ordered ``list`` of ``Region``\ s in components.
    junctions : list
        Ordered ``list`` of ``Junction``\ s in components.

    Methods
    -------
    contains(self, other_sl, junction_tolerance)
        Checks if the given other ``SpliceList`` can be contained within this one.
    """

    def __init__(self, identifier, exons,
                set_left_junction=False, set_right_junction=False):
        r"""
        Parameters
        ----------
        identifier : str
            Accession number if this ``SpliceList`` represents a splice variant,
            otherwise 'pre-mRNA' for precursor mRNAs.
        exons : sequence
            Sequence iterable containing two-element tuples representing exons.
            The two elements are the start and stop positions of this exon,
            respectively. The positions are 1-based and inclusive on both sides.
            The start position should be less than or equal to the stop.
        set_left_junction : bool
            If the 'leftmost' (lower numerical value) nucleotide of the leftmost
            exon should be identified with a ``Junction``.
        set_right_junction : bool
            If the 'rightmost' (higher numerical value) nucleotide of the
            rightmost exon should be identified with a ``Junction``.
        """

        self.identifier = identifier
        exons = sorted(exons, key=lambda pair: pair[0])

        # make sure that none of the exons overlap
        for i in range(0, len(exons)-1):
            this_exon, next_exon = exons[i:i+2]
            if this_exon[1] >= next_exon[0]:
                raise OverlappingExonsError('in {}: {} overlaps {}'.format(
                        identifier, this_exon, next_exon))

        components = []
        for exon in exons:
            components.append(Junction(exon[0], Junction.TYPE_START))
            components.append(Region(*exon))
            components.append(Junction(exon[1], Junction.TYPE_END))
        self.components = components

        if set_left_junction:
            self.components[0].has_complement = True
        else:
            self.components = self.components[1:]
        if set_right_junction:
            self.components[-1].has_complement = True
        else:
            self.components = self.components[:-1]

        self.regions = list(filter(
                lambda x: type(x) == Region, self.components))
        self.junctions = list(filter(
                lambda x: type(x) == Junction, self.components))

    def contains(self, other_sl, junction_tolerance):
        r"""
        Checks if the given SpliceList (other_sl) is contained within this
        SpliceList. In particular,

        1. All regions of other_sl must be fully contained within the regions of
           this SplceList, with tolerance given by junction_tolerance.
        2. All junctions of the other_sl must also be junctions of this
           SplceList, with tolerance given by junction_tolerance. The junctions
           must also match **in order**, i.e. no skips are allowed. This
           identifies splice variants with exon skipping more sensitively.
        """
        # Find the intersect of other_sl.regions and self.regions
        other_regions = set.union(
                *[ set(range(r.start, r.stop+1)) for r in other_sl.regions ])
        this_regions = set.union(
                *[ set(range(r.start, r.stop+1)) for r in self.regions ])
        # Set difference is not associative! For overhangs, take others - this
        overhangs = list(other_regions.difference(this_regions))
        # Check if any other the overhangs > junction_tolerance
        max_overhang_distance = 0
        last_position = -1  # an init value, positions >= 1
        for pos in overhangs:
            # if adjacent to the last position
            if pos - last_position == 1:
                max_overhang_distance += 1
            else:
                max_overhang_distance = 1
            last_position = pos

            if max_overhang_distance > junction_tolerance:
                return False

        # For all junctions in other_sl, check that there is a match here
        match_order = []
        for j in other_sl.junctions:
            has_match = False
            check_range = set(range(
                    j.position-junction_tolerance,
                    j.position+junction_tolerance+1))
            for i, this_j in enumerate(self.junctions):
                if this_j.type == j.type and this_j.position in check_range:
                    if j.complement != None:
                        match_order.append(i)  # preparing to check complements
                    has_match = True
                    continue
            if not has_match:
                return False

        logger.debug(
                '{} contains {}?\n\tself.junctions: {}\n\t'
                'other_sl.junctions: {}'.format(
                    self.identifier, other_sl.identifier,
                    [str(j) for j in self.junctions],
                    [str(j) for j in other_sl.junctions]))
        # Check if the transcript's complementing Junctions are also
        # complementing in the annotation.
        for i in range(0, len(match_order)-1, 2):
            logger.debug('{} contains {}? match_order[{}:{}] == [{}, {}]'
                    .format(self.identifier, other_sl.identifier, i, i+2,
                    match_order[i], match_order[i+1]))
            if abs(match_order[i] - match_order[i+1]) > 1:
                return False

        return True

    def __str__(self):
        return '{} <{}>'.format(
                self.identifier,
                ', '.join(map(str, self.components)))

    @staticmethod
    def join(sl1, sl2):
        r"""
        Combines two ``SpliceList``\ s. The ``Region``\ s of the new
        ``SpliceList`` are the union of pre-existing ``Region``\ s. The
        ``Junction`` s are the union of pre-existing ``Juntion``\ s, with
        ``Junction`` which lie within an unbroken ``Region`` removed.

        Attributes
        ----------
        sl1 : SpliceList
        sl2 : SpliceList

        Returns
        -------
        SpliceList
        """
        identifier = sl1.identifier if sl1.identifier == sl2.identifier else \
                '{},{}'.format(sl1.identifier, sl2.identifier)
        mapped = [set(range(r.start, r.stop+1)) \
                for r in sl1.regions + sl2.regions]
        mapped = sorted(list(set.union(*mapped)))

        exons = []
        start = None
        last = 0
        for pos in mapped:
            if start == None:
                start = pos
                last = pos
            elif pos - last == 1:
                last += 1
            else:
                exons.append([start, last])
                start = pos
                last = pos
        else:
            exons.append([start, pos])

        regions = [Region(*coords) for coords in exons]
        _junctions = set(sl1.junctions).union(sl2.junctions)
        junctions = []
        for j in _junctions:
            if j.type == Junction.TYPE_START and not j.position-1 in mapped:
                junctions.append(j)
            elif j.type == Junction.TYPE_END and not j.position+1 in mapped:
                junctions.append(j)

        components = regions + junctions
        components = sorted(components,
                key=lambda c: c.start if type(c) == Region else c.position)
        regions = list(filter(lambda c: type(c) == Region, components))
        junctions = list(filter(lambda c: type(c) == Junction, components))

        new_sl = SpliceList('', [])
        new_sl.identifier = identifier
        new_sl.components = components
        new_sl.junctions = junctions
        new_sl.regions = regions

        return new_sl

    @staticmethod
    def join_many(*sls):
        r"""
        Combines two or more SpliceLists by calling the ``SpliceList.join``
        static method.

        Attributes
        ----------
        sls : iterable
            Sequence iterable of more than one ``SpliceList``\ s.

        Returns
        -------
        SpliceList
        """
        joined = sls[0]
        for sl in sls[1:]:
            joined = SpliceList.join(joined, sl)

        logger.debug('Joined transcripts to\n\t{}'.format(joined))
        return joined


class Region:
    r"""
    A nucleotide region.

    Attributes
    ----------
    start : int
        1-based and inclusive start position of this region. Should be less than
        or equal to the ``stop``.
    stop : int
        1-based and inclusive stop position of this region. Should be greater
        than or equal to the ``start``.
    """

    def __init__(self, start, stop):
        r"""
        Parameters
        ----------
        start : int
            1-based and inclusive start position of this region. Should be less
            than or equal to the ``stop``.
        stop : int
            1-based and inclusive stop position of this region. Should be
            greater than or equal to the ``start``.
        """
        self.start = start
        self.stop = stop

    def __str__(self):
        return 'region {}-{}'.format(self.start, self.stop)


class Junction:
    r"""
    A splice junction.

    Attributes
    ----------
    position : int
        1-based position of the exon of this junction.
    type : int
        If this junction corresponds to the start or end of an exon, based on a
        numerically increasing order of nucleotide position (**not** concerned
        with 3' or 5' directionality).
    TYPE_START : int
        Constant, class attribute for ``type`` representing a junction type which
        leads into an exon, based on a numerically increasing order of
        nucleotide position.
    TYPE_END : int
        Constant, class attribute for ``type`` representing a junction type which
        leads into an intron, based on a numerically increasing order of
        nucleotide position.
    """

    TYPE_START = 0
    TYPE_END = 1

    def __init__(self, position, _type):
        r"""
        Parameters
        ----------
        position : int
            1-based position of the exon of this junction.
        type : int
            If this junction corresponds to the start or end of an exon, based on a
            numerically increasing order of nucleotide position (**not** concerned
            with 3' or 5' directionality).
        complement : Junction
            The ``Junction`` of complement type in the read group.  E.g. if this
            `Junction` is `TYPE_END`, then its ``complement`` is the
            ``TYPE_START`` ``Junction`` immediately adjacent to it on the mature
            mRNA read.
        """
        self.position = position
        self.type = _type

        self.complement = None
        self.has_complement = False

    def __str__(self):
        return 'junction ({}) at {}{}'.format(
                'start' if self.type == self.TYPE_START else 'end',
                self.position,
                ' with complement at {}'.format(self.complement.position)
                    if self.complement != None else '')

class OverlappingExonsError(Exception):
    r"""
    ``Exception`` raised if two regions in a ``SpliceList`` are found to overlap
    each other.
    """
    pass
