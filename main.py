#!/usr/bin/env python3
import argparse
import logging
import random
from collections import Counter

logger = logging.getLogger(__name__)


def _fastq_parser(path):
    header = seq = ''

    with open(path) as f:
        for l in f:
            if l.startswith('>'):
                if seq:
                    yield header, seq.upper()
                    seq = ''
                header = l.lstrip('>').strip()
            else:
                seq += l.rstrip()

    if header:
        yield header, seq.upper()


def _main():
    parser = argparse.ArgumentParser()

    parser.add_argument("fasta", nargs=1)
    parser.add_argument("-a", "--alpha", type=float, default=0.)
    parser.add_argument("--ignore-n-bases", action="store_false")
    parser.add_argument("--seed", nargs='?')
    parser.add_argument("-n", "--peaks", type=int, default=20000)
    parser.add_argument("-w", "--width", type=int, default=100)
    parser.add_argument("-r", "--reads", type=int, default=10000000)
    parser.add_argument("-d", "--distance", type=int, default=200)
    outfmt = parser.add_mutually_exclusive_group(required=True)
    outfmt.add_argument("-q", "--fastq", action="store_true")
    outfmt.add_argument("-s", "--sam", action="store_true")

    args = parser.parse_args()

    if args.seed:
        random.seed(args.seed)

    name2seq = {}
    name2len = {}
    for header, seq in _fastq_parser(args.fasta[0]):
        name2seq[header] = seq

        c = Counter(seq)
        if args.ignore_n_bases:
            del c['N']
        name2len[header] = sum(c.values())

        logger.info("Load '{}', len = {:,}".format(header, name2len[header]) + ' ' + repr(c))


if __name__ == "__main__":
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)
    _main()
