#!/usr/bin/env python3
import argparse
import logging
import random
import fnmatch
from collections import Counter
from itertools import groupby, chain
import sys
import functools

from intervaltree import Interval, IntervalTree

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


class Chromosome:
    COUNT_N = str.maketrans("ATGC", 'Y' * 4)

    def __init__(self, seq, fraglen, ignore_n=True):
        self.seq = seq
        self.chromlen = self.len = len(seq)
        self.peak_regions = IntervalTree()
        self.blacklist = IntervalTree()
        self.fraglen = fraglen
        if ignore_n:
            pos = 0
            for k, g in groupby(self.seq.translate(self.COUNT_N)):
                l = sum(1 for _ in g)
                if k == 'N':
                    self.blacklist.add(Interval(pos, pos + l))
                    self.len -= l
                pos += l

    def choose_peak_regions(self, n):
        while len(self.peak_regions) < n:
            pos = random.randrange(self.chromlen - self.fraglen)
            peak = Interval(pos, pos + self.fraglen)
            if not self.blacklist.overlap(peak) and not self.peak_regions.overlap(peak):
                self.peak_regions.add(peak)

    def _get_read_from_fragment(self, frag, width, readlen):
        positive_strand = bool(random.getrandbits(1))
        if positive_strand:
            pos = random.randrange(frag.begin, frag.begin + width)
        else:
            pos = random.randrange(frag.end - width, frag.end)
            pos -= readlen - 1

        return pos, positive_strand

    def get_reads_from_peaks(self, width, readlen, n):
        peaks = tuple(self.peak_regions)
        if not peaks:
            raise StopIteration

        for _ in range(n):
            peak = random.choice(peaks)
            yield self._get_read_from_fragment(peak, width, readlen)

    def get_reads_as_background(self, width, readlen, n):
        for _ in range(n):
            pos = random.randrange(0, self.chromlen - self.fraglen)
            fragment = Interval(pos, pos + self.fraglen)
            if not self.blacklist.overlap(fragment):
                yield self._get_read_from_fragment(fragment, width, readlen)


class Formatter:
    COMPLEMENT = str.maketrans("ATGC", "TACG")

    @staticmethod
    def _make_seq_name(chrname, pos, is_pos_strand, is_peak):
        return "{}_{}_{}_{}".format(
            chrname,
            pos,
            '+' if is_pos_strand else '-',
            "peak" if is_peak else "background"
        )

    @staticmethod
    def revcom(func):
        @functools.wraps(func)
        def _revcomp(seq, is_pos_strand, *args, **kwargs):
            if not is_pos_strand:
                seq = seq.translate(Formatter.COMPLEMENT)[::-1]
            return func(seq, is_pos_strand, *args, **kwargs)
        return _revcomp


class FASTQFormatter(Formatter):
    @Formatter.revcom
    def format(self, chrname, pos, is_pos_strand, seq, is_peak):
        print('@' + self._make_seq_name(chrname, pos, is_pos_strand, is_peak),
              seq, '+' if is_pos_strand else '-', 'I' * len(seq), sep='\n')


class SAMFormatter(Formatter):
    MAPQ = 42

    def __init__(self, names, lengths):
        print("@HD", "VN:1.0", "SO:coordinate", sep='\t')
        for n, l in zip(names, lengths):
            print("@SQ", "SN:" + n, "LN:{}".format(l), sep='\t')
        print("@PG", "ID:Python{0}.{1}.{2}".format(*sys.version_info), "PN:" + __name__, "CL:" + ' '.join(sys.argv))

    @Formatter.revcom
    def format(self, chrname, pos, is_pos_strand, seq, is_peak):
        print(
            self._make_seq_name(chrname, pos, is_pos_strand, is_peak),
            0 if is_pos_strand else 16,
            chrname,
            pos + 1,
            self.MAPQ,
            "{}M".format(len(seq)),
            '*',
            0,
            0,
            seq,
            'I' * len(seq),
            sep='\t'
        )


def _main():
    parser = argparse.ArgumentParser()

    parser.add_argument("fasta", nargs=1)
    parser.add_argument("-a", "--alpha", type=float, default=0.01)
    parser.add_argument("--ignore-n-bases", action="store_false")
    parser.add_argument("-i", "--include-chr", nargs='+')
    parser.add_argument("--seed", nargs='?')
    parser.add_argument("-n", "--peaks", type=int, default=20000)
    parser.add_argument("-w", "--width", type=int, default=100)
    parser.add_argument("-r", "--reads", type=int, default=10000000)
    parser.add_argument("-d", "--distance", type=int, default=200)
    parser.add_argument("-l", "--readlen", type=int, default=50)
    outfmt = parser.add_mutually_exclusive_group(required=True)
    outfmt.add_argument("-q", "--fastq", action="store_true")
    outfmt.add_argument("-s", "--sam", action="store_true")
    parser.add_argument("-v", "--verbose", action="store_true")

    args = parser.parse_args()
    fraglen = args.width * 2 + args.distance

    if args.seed:
        random.seed(args.seed)
    if args.verbose:
        logger.setLevel(logging.INFO)

    chrnames = []
    chrlengths = []
    name2chr = {}
    for header, seq in _fastq_parser(args.fasta[0]):
        if args.include_chr and not any(fnmatch.fnmatch(header, pattern) for pattern in args.include_chr):
            continue
        chrnames.append(header)
        name2chr[header] = c = Chromosome(seq, fraglen, args.ignore_n_bases)
        chrlengths.append(c.chromlen)
        logger.info("Load '{}', len = {:,}".format(header, c.chromlen))

    #
    genomelen = sum(c.chromlen for c in name2chr.values())
    if fraglen * args.peaks > genomelen / 2:
        logger.error("Too many peak regions.")
        sys.exit(1)
    peak_sizes = Counter(random.choices(chrnames, weights=chrlengths, k=args.peaks))
    logger.info("generate peaks: {}".format(peak_sizes))
    for name in chrnames:
        name2chr[name].choose_peak_regions(peak_sizes[name])

    #
    if args.fastq:
        formatter = FASTQFormatter()
    elif args.sam:
        formatter = SAMFormatter(chrnames, chrlengths)
    else:
        raise AssertionError

    #
    sample_sizes = Counter(random.choices(chrnames, weights=chrlengths, k=args.reads))
    logger.info("generate reads: {}".format(sample_sizes))
    for name in chrnames:
        chrom = name2chr[name]
        reads_in_peaks = int(sample_sizes[name] * args.alpha)
        reads_in_bg = sample_sizes[name] - reads_in_peaks
        reads = sorted(
            chain(
                ((pos, strand, False) for pos, strand in chrom.get_reads_as_background(args.width, args.readlen, reads_in_bg)),
                ((pos, strand, True) for pos, strand in chrom.get_reads_from_peaks(args.width, args.readlen, reads_in_peaks))
            ),
            key=lambda x: x[0]
        )
        for pos, is_pos_strand, is_peak in reads:
            seq = chrom.seq[pos:(pos + args.readlen)]
            formatter.format(name, pos, is_pos_strand, seq, is_peak)


if __name__ == "__main__":
    logger.addHandler(logging.StreamHandler())
    _main()
