import logging

import cli
from transcript import get_transcripts
from annotations import Annotations

VERSION = '0.0.1-SNAPSHOT'

def main():
    r"""
    Runs trcls.
    """
    parser = cli.get_parser()
    args = parser.parse_args()
    logger = setup_logging(args)

    if args.version:
        print('trcls {}'.format(VERSION))
        exit(0)

    if args.alignment == None or args.features == None:
        logger.error(
                'Both SAM alignment and GTF annotation files must be provided')
        parser.print_help()
        exit(1)

    with open(args.features) as features_file:
        annotations = Annotations(features_file)

    with open(args.alignment) as alignment_file:
        alignments = alignment_file.readlines()

    headers = filter(lambda l: l.startswith('@'), alignments)
    headers = map(str.strip, headers)
    alignments = filter(lambda l: not l.startswith('@'), alignments)

    transcripts = get_transcripts(
            alignments, args.skip_tolerance, args.map_tolerance)

    print('\n'.join(headers))
    for transcript in transcripts:
        transcript.annotate(annotations, args.junction_tolerance)
        print(transcript)

def setup_logging(args):
    r"""
    Initiates the logging library.

    Parameters
    ----------
    args : Namespace
        `Namespace` object from `ArgumentParser.parse_args`

    Returns
    -------
    `logging.Logger` object.
    """
    logging.basicConfig(format='[%(asctime)s/%(levelname)s] %(message)s')
    logger = logging.getLogger('trcls')

    if args.quiet:
        logger.setLevel(logging.ERROR)
    elif args.verbose:
        logger.setLevel(logging.INFO)
    elif args.very_verbose:
        logger.setLevel(logging.DEBUG)

    return logger

if __name__ == '__main__':
    main()
