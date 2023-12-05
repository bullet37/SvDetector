from utils import *
from utilsPE import *
import pysam
import argparse
import importlib
import os
def detect():
    file_url = r" "
    c = Config()
    c.files = [file_url]
    svs = []

    # Open BAM file using pysam
    samfile = pysam.AlignmentFile(c.files[0], "rb")
    chromosome_number = len(samfile.references)
    c.nchr = chromosome_number

    # Get the header
    hdr = samfile.header
    # print("Chromosome number: ", len(hdr['SQ'])) # len(hdr['SQ'])== len(samfile.references)
    validRegions = [IntervalSet() for _ in range(len(samfile.references))]
    if not _parseExcludeIntervals(c, samfile, validRegions):
        print("Delly couldn't parse exclude intervals!")
        exit(1)

    # Create library objects
    sampleLib = [LibraryInfo() for _ in c.files]
    getLibraryParams(c, validRegions, sampleLib)

    if not c.hasVcfFile:
        srSVs = []
        srStore = [{} for _ in range(len(hdr['SQ']))]
        bam,svs= scanPE(c, validRegions, svs, srSVs, srStore, sampleLib)

    for i, item in enumerate(svs):
      item.print_info()

def main():
    parser = argparse.ArgumentParser(description="SvDetector")

    parser.add_argument('--ref_path', type=str, default="./", help='Reference fatsq path')
    parser.add_argument('--bamfile_path', type=str, required=True, help='Bamfile path')
    parser.add_argument('--output_type', type=str, required=True, help='Output vcf type')
    # parser.add_argument('-- ', type=float, default=0.001, help='Learning Rate')

    args = parser.parse_args()

    ref_path = args.ref_path
    bamfile_path = args.bamfile_path
    output_type = args.output_type

    detect()

if __name__ == "__main__":
    main()