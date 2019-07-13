import argparse

DEFAULT_JUNCTION = 20
DEFAULT_MAP = 10
DEFAULT_SKIP = 20

def get_parser():
    r"""
    Parse the command line arguments

    Returns
    -------
    Namespace
        Arguments and options passed to the command line for the execution of
        this programe.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('alignment', nargs='?', help='a SAM alignment file')
    parser.add_argument('features', nargs='?', help='a GTF annotation file')

    parser.add_argument('--map-tolerance', '-m',
            type=int, default=DEFAULT_MAP,
            help=(
                    'number of nt which must map to reference in sequence to '
                    'be considered a feature, e.g. exon (default: {})'.format(
                            DEFAULT_MAP)))
    parser.add_argument('--skip-tolerance', '-s',
            type=int, default=DEFAULT_SKIP,
            help=(
                    'number of non-splice reference skips allowed, e.g. '
                    ' in deletions (default: {})'.format(DEFAULT_SKIP)))
    parser.add_argument('--junction-tolerance', '-o',
            type=int, default=DEFAULT_JUNCTION,
            help=(
                    'number of nt exon-intron junctions allowed to mismatch by '
                    '(default: {})'.format(DEFAULT_JUNCTION)))

    parser.add_argument('--quiet', '-q', action='store_true',
            help='show only error (ERROR) messages')
    parser.add_argument('--verbose', '-v', action='store_true',
            help='show verbose (INFO, WARNING, ERROR) messages')
    parser.add_argument('--very-verbose', '-vv', action='store_true',
            help='show very verbose (DEBUG, INFO, WARNING, ERROR) messages')
    parser.add_argument('--version', action='store_true',
            help='Print the version number')

    return parser
