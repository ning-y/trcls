import logging
from splice_list import SpliceList

logger = logging.getLogger('trcls')

class Annotations:
    r"""
    A set of known splice variants and precursor mRNA.

    Attributes
    ----------
    variants : list
        ``list`` of ``SpliceList``\ s representing known splice variants, and
        precursor mRNA.

    Methods
    -------
    get_annotations(transcript, junction_tolerance)
        Given a ``Transcript``, return a ``list`` of ``str`` identifiers of splice
        variants or precursor mRNA which the transcript could have originated
        from.
    """

    def __init__(self, gtf_file):
        r"""
        Parameters
        ----------
        gtf_file : file
            Annotation file describing all splice variants found in the region
            of the alignment.
        """
        # lines is a list of exons only
        lines = gtf_file.readlines()
        lines = map(lambda line: line.split('\t'), lines)
        lines = filter(lambda line: line[2] == 'exon', lines)
        lines = list(lines)

        lines2 = []  # contains tuple<start, stop, transcript_id>
        for line in lines:
            start, stop = int(line[3]), int(line[4])
            transcript_id = line[-1].split(';')
            transcript_id = next(filter(
                    lambda field: 'transcript_id' in field, transcript_id))
            transcript_id = transcript_id.replace('transcript_id', '')
            transcript_id = transcript_id.replace('"', '')
            transcript_id = transcript_id.strip()
            lines2.append((start, stop, transcript_id))

        precursor_start = None
        precursor_end = None
        variants = []
        exons = []
        current_variant = lines2[0][2]
        for triple in lines2:
            # Tracking exons of the same splice variant
            if triple[-1] == current_variant:
                exons.append(triple[0:2])
            # A new splice variant has appeared
            else:
                # Add it to variants as a SpliceList
                should_reverse = exons[0][0] > exons[0][1]
                if should_reverse:
                    exons = list(map(sorted, exons))
                    exons = sorted(exons, key=lambda pair: pair[0])
                variants.append(SpliceList( current_variant, exons,
                        set_left_junction=True, set_right_junction=True))

                precursor_start = exons[0][0] if precursor_start == None \
                        else min(precursor_start, exons[0][0])
                precursor_end = exons[-1][1] if precursor_end == None \
                        else max(precursor_end, exons[-1][1])

                # Reset
                exons = [triple[0:2]]
                current_variant = triple[2]
        else:
            # Clean-up the last splice variant
            should_reverse = exons[0][0] > exons[0][1]
            if should_reverse:
                exons = list(map(sorted, exons))
                exons = sorted(exons, key=lambda pair: pair[0])
            variants.append(SpliceList(current_variant, exons,
                    set_left_junction=True, set_right_junction=True))

            precursor_start = exons[0][0] if precursor_start == None \
                    else min(precursor_start, exons[0][0])
            precursor_end = exons[-1][1] if precursor_end == None \
                    else max(precursor_end, exons[-1][1])
            variants.append(
                    SpliceList('pre-mRNA', [(precursor_start, precursor_end)],
                            set_left_junction=True, set_right_junction=True))

        self.variants = variants
        logger.debug('Instantiated annotations with\n{}'.format(
                '\n\t'.join([str(sl) for sl in self.variants])))

    def get_annotations(self, transcript, junction_tolerance):
        r"""
        Find the possible origins of the transcript from this annotation.

        Parameters
        ----------
        transcript : Transcript
            The query transcript.
        junction_tolerance : int
            Number of nucleotides which the transcript junctions are allowed to
            differ in position from a splice variant's junction such that the
            transcript can still be considered as being from the splice variant.
            This is particularly useful for transcripts aligned with non-splice
            aware aligners that have a tendency to slightly misalign transcripts
            near the splice junction.

        Returns
        -------
        list
            ``list`` of ``str`` identifiers of the possible splice variants or
            precursor mRNA which this transcript could have originated from.
        """
        # Get list of all splice_lists in transcript
        transcript_splice_lists = list(map(
                lambda t: t.splice_list, transcript.segments))
        # Merge the splice_lists
        transcript_splice_lists_merged = SpliceList.join_many(
                *transcript_splice_lists)
        # Get True/False for possible matches in self.variants
        matches = map(lambda v:
                v.contains(transcript_splice_lists_merged, junction_tolerance),
                self.variants)
        possible_annotations = []
        for variant, contains_transcript in zip(self.variants, matches):
            if contains_transcript:
                possible_annotations.append(variant.identifier)

        return possible_annotations
