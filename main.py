from scannerPE import *
from interfaceUtils import *
import pysam
import argparse
import os
import random
from VCFwriter import *


def parse_arguments():
    parser = argparse.ArgumentParser(description="SvDetector command line options")
    # Generic options
    parser.add_argument("-t", "--svtype", default="ALL", help="SV type to compute [DEL, INS, DUP, INV, BND, ALL]")
    parser.add_argument("-g", "--genome", help="genome fasta file")
    parser.add_argument("-x", "--exclude", help="file with regions to exclude")
    parser.add_argument("-o", "--outfile", help="BCF output file")

    parser.add_argument("-q", "--map-qual", type=int, default=1, help="min paired-end (PE) mapping quality")
    parser.add_argument("-r", "--qual-tra", type=int, default=20, help="min PE quality for translocation")
    parser.add_argument("-s", "--mad-cutoff", type=int, default=9,
                        help="insert size cutoff, median+s*MAD (deletions only)")
    parser.add_argument("-c", "--minclip", type=int, default=25, help="min clipping length")
    parser.add_argument("-z", "--min-clique-size", type=int, default=2, help="min. PE/SR clique size")
    parser.add_argument("-m", "--minrefsep", type=int, default=25, help="min reference separation")
    parser.add_argument("-n", "--maxreadsep", type=int, default=40, help="max read separation")

    parser.add_argument("-v", "--vcffile", help="input VCF/BCF file for genotyping")
    parser.add_argument("-u", "--geno-qual", type=int, default=5, help="min mapping quality for genotyping")
    parser.add_argument("-d", "--dump", help="gzipped output file for SV-reads (optional)")
    parser.add_argument("-a", "--max-geno-count", type=int, default=250,
                        help="max. number of reads aligned for SR genotyping")

    # Hidden options
    parser.add_argument("-j", "--pruning", type=int, default=1000, help="PE graph pruning cutoff")
    parser.add_argument("-w", "--cons-window", type=int, default=100, help="consensus window")
    parser.add_argument("input_files", nargs='*', help="input files")

    args = parser.parse_args()
    if not args.input_files or not args.genome:
        print("\nUsage: SvDetector [OPTIONS] -g <ref.fa> <sample1.sort.bam> <sample2.sort.bam> ...")
        parser.print_help()
        exit(1)

    return args


def main():
    args = parse_arguments()
    config = Config()
    config.setup(args)
    svs = []

    config.outfile = "test.vcf"
    # Open BAM file using pysam
    samfile = pysam.AlignmentFile(config.files[0], "rb")
    config.nchr = len(samfile.references)

    # Get the header
    hdr = samfile.header
    # print("Chromosome number: ", len(hdr['SQ'])) # len(hdr['SQ'])== len(samfile.references)
    validRegions = [IntervalSet() for _ in range(len(samfile.references))]
    if not parseExcludeIntervals(config, samfile, validRegions):
        print("SvDetector couldn't parse exclude intervals!")
        exit(1)

    # Create library reference
    sampleLib = [LibraryInfo() for _ in config.files]
    get_library_params(config, validRegions, sampleLib)

    if not config.hasVcfFile:
        srSVs = []
        srStore = [{} for _ in range(len(hdr['SQ']))]
        scanPE(config, validRegions, svs, srSVs, srStore, sampleLib)

    for i, item in enumerate(svs):
        item.print_info()

    svs.sort(key=cmp_to_key(sort_svs))

    # todo: global Annotation
    # Fake data:
    class JunctionCount:
        def __init__(self):
            self.ref = [random.randint(1, 10) for _ in range(5)]
            self.alt = [random.randint(1, 10) for _ in range(5)]

    class ReadCount:
        def __init__(self, l=0, m=0, r=0):
            self.leftRC = l
            self.rc = m
            self.rightRC = r

    class SpanningCount:
        def __init__(self):
            self.ref = [random.randint(1, 10) for _ in range(5)]
            self.alt = [random.randint(1, 10) for _ in range(5)]

    jctMap = [[JunctionCount() for _ in range(2)] for _ in range(2)]
    read_count_array = [
        [ReadCount(random.randint(0, 10), random.randint(0, 10), random.randint(0, 10)) for _ in range(2)] for _ in
        range(2)]
    spanning_count_array = [[SpanningCount() for _ in range(2)] for _ in range(2)]

    vcf_output(config, svs, jctMap, read_count_array, spanning_count_array)


if __name__ == "__main__":
    main()
