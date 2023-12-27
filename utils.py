# Config
from typing import List, Dict, Set
from functools import cmp_to_key
import numpy as np
import pysam
import sys, os

SVT_TRANS = 5


class Config:
    def __init__(self):
        self.min_map_qual = 1
        self.min_geno_qual = 5
        self.nchr = 0
        self.minimum_flank_size = 13
        self.indel_size = 1000
        self.min_cons_window = 100
        self.graph_pruning = 1000
        self.min_ref_sep = 25
        self.max_read_sep = 40
        self.min_clip = 25
        self.max_geno_read_count = 250
        self.min_clique_size = 2
        self.flank_quality = 0.95

        self.has_dump_file = False
        self.has_exclude_file = False
        self.has_vcf_file = False
        self.svtset = set()
        self.aliscore = DnaScore(5, -4, -10, -1)
        self.dumpfile = ''
        self.outfile = ''
        self.vcffile = ''
        self.genome = ''
        self.exclude = ''
        self.files = []
        self.sample_name = []

        self.min_tra_qual = 20
        self.mad_cutoff = 9
        self.mad_normal_cutoff = 5

    def print_info(self):
        print("Config Information:")
        for key, value in self.__dict__.items():
            print(f"{key}: {value}")

    def _set_svtset(self, svtype):
        SVT_TRANS = 5
        if svtype == "ALL":
            return True
        else:
            svtypes = svtype.split(',')
            for tok in svtypes:
                if tok == "DEL":
                    self.svtset.add(2)
                elif tok == "INS":
                    self.svtset.add(4)
                elif tok == "DUP":
                    self.svtset.add(3)
                elif tok == "INV":
                    self.svtset.add(0)
                    self.svtset.add(1)
                elif tok == "INV_3to3":
                    self.svtset.add(0)
                elif tok == "INV_5to5":
                    self.svtset.add(1)
                elif tok == "BND":
                    self.svtset.update(range(SVT_TRANS, SVT_TRANS + 4))
                elif tok == "BND_3to3":
                    self.svtset.add(SVT_TRANS + 0)
                elif tok == "BND_5to5":
                    self.svtset.add(SVT_TRANS + 1)
                elif tok == "BND_3to5":
                    self.svtset.add(SVT_TRANS + 2)
                elif tok == "BND_5to3":
                    self.svtset.add(SVT_TRANS + 3)
                else:
                    return False
            return True

    def check(self):
        if self.exclude:
            if not (os.path.exists(self.exclude) and os.path.isfile(self.exclude) and os.path.getsize(
                    self.exclude) > 0):
                print(f"Exclude file is missing: {self.exclude}", file=sys.stderr)
                sys.exit(1)
            self.hasExcludeFile = True

        if not (os.path.exists(self.genome) and os.path.isfile(self.genome) and os.path.getsize(self.genome) > 0):
            print(f"Reference file is missing: {self.genome}", file=sys.stderr)
            sys.exit(1)

        if self.vcffile:
            if not (os.path.exists(self.vcffile) and os.path.isfile(self.vcffile) and os.path.getsize(
                    self.vcffile) > 0):
                print(f"Input VCF/BCF file is missing: {self.vcffile}", file=sys.stderr)
                sys.exit(1)
            try:
                with pysam.VariantFile(self.vcffile, "r") as ifile:
                    hdr = ifile.header
            except Exception as e:
                print(f"Fail to open file {self.vcffile}: {e}", file=sys.stderr)
                sys.exit(1)
            self.hasVcfFile = True

        # if not self.outfile:
        #     self.outfile = "-"
        # else:
        #     if self.outfile != "-":
        #         if not _outfileValid(self.outfile):  # Todo
        #             sys.exit(1)
        if self.min_geno_qual < 5:
            self.min_geno_qual = 5

        for file_c, file_path in enumerate(self.files):
            if not (os.path.exists(file_path) and os.path.isfile(file_path) and os.path.getsize(file_path) > 0):
                print(f"Alignment file is missing: {file_path}", file=sys.stderr)
                sys.exit(1)
            try:
                samfile = pysam.AlignmentFile(file_path, "rb")
                hdr = samfile.header
            except IOError as e:
                print(f"Fail to open file {file_path}: {e}", file=sys.stderr)
                sys.exit(1)
            n_targets = len(hdr.references)
            if self.nchr == 0:
                self.nchr = n_targets
            elif self.nchr != n_targets:
                print("BAM files have different number of chromosomes!", file=sys.stderr)
                sys.exit(1)
            with pysam.FastaFile(self.genome) as fai:
                for ref_name in hdr.references:
                    if ref_name not in fai.references:
                        print(f"BAM file chromosome {ref_name} is NOT present in your reference file {self.genome}",
                              file=sys.stderr)
                        sys.exit(1)
            samfile.close()

    def setup(self, args):
        self.genome = args.genome
        self.exclude = args.exclude
        self.outfile = args.outfile
        self.min_map_qual = args.map_qual
        self.min_tra_qual = args.qual_tra
        self.mad_cutoff = args.mad_cutoff
        self.min_clip = args.min_clip
        self.min_ref_sep = args.min_ref_sep
        self.max_read_sep = args.max_read_sep
        self.vcffile = args.vcffile
        self.min_geno_qual = args.geno_qual
        self.dumpfile = args.dump
        self.max_geno_read_count = args.max_geno_count
        self.files = args.input_files
        self.graph_pruning = args.pruning

        self.aliscore = DnaScore(5, -4, -10, -1)
        self.hasDumpFile = bool(args.dump)
        self.min_clique_size = max(2, args.min_clique_size)
        self.min_tra_qual = max(self.min_map_qual, args.qual_tra)

        if not self._set_svtset(args.svtype):
            print("Please specify a valid SV type, i.e., -t INV or -t DEL,INV without spaces.")
            sys.exit(1)

        for file_c, file_path in enumerate(self.files):
            try:
                samfile = pysam.AlignmentFile(file_path, "rb")
                hdr = samfile.header
            except IOError as e:
                print(f"Fail to open file {file_path}: {e}", file=sys.stderr)
                sys.exit(1)
            sample_name = "unknown"
            sample_name = self._get_SM_tag(hdr.text, os.path.splitext(os.path.basename(file_path))[0])
            self.sample_name.append(sample_name)
            samfile.close()

    def _get_SM_tag(self, header, file_name):
        sm_identifiers = set()
        rg_present = False
        lines = header.split('\n')
        for line in lines:
            if line.startswith("@RG"):
                key_values = line.split("\t")
                for key_value in key_values:
                    if ":" in key_value:
                        field, value = key_value.split(":", 1)
                        if field == "SM":
                            rg_present = True
                            sm_identifiers.add(value)

        if not rg_present:
            return file_name
        elif len(sm_identifiers) == 1:
            return next(iter(sm_identifiers))
        elif len(sm_identifiers) > 1:
            print("Warning: Multiple sample names (@RG:SM) present in the BAM file!", file=sys.stderr)
            return next(iter(sm_identifiers))


class DnaScore:
    def __init__(self, match=5, mismatch=-4, gap_open=-10, gap_extension=-1):
        self.match = match
        self.mismatch = mismatch
        self.gap_open = gap_open
        self.gap_extension = gap_extension
        self.inf = 1000000

    def print_info(self):
        print("DnaScore Information:")
        for key, value in self.__dict__.items():
            print(f"{key}: {value}")


class LibraryInfo:
    def __init__(self, read_size=0, median=0, mad=0, min_normal_isize=0,
                 min_isize_cutoff=0, max_normal_isize=0, max_isize_cutoff=0, abnormal_pairs=0):
        self.read_size = read_size  # median read size
        self.median = median
        self.mad = 0
        self.min_normal_isize = 0
        self.min_isize_cutoff = 0
        self.max_normal_isize = 0
        self.max_isize_cutoff = 0
        self.abnormal_pairs = 0

    def print_info(self):
        print("Library Info:")
        for key, value in self.__dict__.items():
            print(f"{key}: {value}")


def get_library_params(c: Config, validRegions, sampleLib):
    samfiles = [pysam.AlignmentFile(f, 'r') for f in c.files]
    for f in c.files:
        pysam.index(f)
    for file_c, file_name in enumerate(c.files):
        max_alignments_screened = 10000000
        max_alignments_num = 1000000
        min_alignments_num = 1000
        alignment_count = 0
        processed_pairs_num = 0
        processed_reads_num = 0
        rplus = 0
        nonrplus = 0
        vecISize = []
        readSize = []

        libCharacterized = False
        samfile = samfiles[file_c]
        for refIndex, refname in enumerate(samfile.references):
            if not validRegions[refIndex]:
                continue
            for vRIt in validRegions[refIndex]:
                for rec in samfile.fetch(refname, vRIt[0], vRIt[1]):
                    # 是否为次要的对齐 rec.is_secondary， rec.tlen = insert size
                    if not rec.is_read2 and rec.infer_read_length() < 65000:
                        if rec.is_secondary or rec.is_qcfail or rec.is_duplicate or rec.is_supplementary or rec.is_unmapped:
                            continue
                        if alignment_count > max_alignments_screened or (
                                processed_reads_num >= max_alignments_num and processed_pairs_num == 0) or processed_pairs_num >= max_alignments_num:
                            libCharacterized = True
                            break
                        alignment_count += 1

                        if processed_reads_num < max_alignments_num:
                            readSize.append(rec.infer_read_length())
                            processed_reads_num += 1

                        if rec.is_paired and not rec.mate_is_unmapped and rec.rname == rec.mrnm:
                            if processed_pairs_num < max_alignments_num:
                                vecISize.append(abs(rec.tlen))
                                if get_sv_type(rec) == 2:
                                    rplus += 1
                                else:
                                    nonrplus += 1
                                processed_pairs_num += 1
                if libCharacterized:
                    break
            if libCharacterized:
                break

        if processed_reads_num >= min_alignments_num:
            readSize.sort()
            sampleLib[file_c].read_size = readSize[len(readSize) // 2]  # Median

        if processed_pairs_num >= min_alignments_num:
            vecISize.sort()
            median = vecISize[len(vecISize) // 2]
            absDev = [abs(i - median) for i in vecISize]
            absDev.sort()
            mad = absDev[len(absDev) // 2]

            if 50 <= median <= 100000:
                if rplus < nonrplus:
                    print(f"Warning: Sample has a non-default paired-end layout! File: {file_name}")
                    print(
                        "The expected paired-end orientation is ---Read1---> <---Read2--- which is the default illumina paired-end layout.")
                else:
                    lib = sampleLib[file_c]
                    lib.median = median
                    lib.mad = mad  # Median Absolute Deviation
                    lib.max_normal_isize = median + (c.mad_normal_cutoff * mad)  # 正常范围threshold
                    lib.min_normal_isize = max(median - (c.mad_normal_cutoff * mad), 0)
                    lib.max_isize_cutoff = max(median + (c.mad_cutoff * mad),
                                               max(2 * lib.read_size, 500))  # 一个异常插入片段的阈值
                    lib.min_isize_cutoff = max(median - (c.mad_cutoff * mad), 0)


class StructuralVariantRecord:
    def __init__(self, chr=0, sv_start=0, chr2=0, sv_end=0, ciposlow=0, ciposhigh=0, ciendlow=0, ciendhigh=0,
                 sr_support=0, pe_support=0, sr_map_quality=0, mapq=0, insLen=0, svt=-1, id=0):
        self.chr = chr
        self.sv_start = sv_start
        self.chr2 = chr2
        self.sv_end = sv_end
        self.ciposlow = ciposlow
        self.ciposhigh = ciposhigh
        self.ciendlow = ciendlow
        self.ciendhigh = ciendhigh

        self.sr_support = sr_support
        self.sr_map_quality = sr_map_quality
        self.pe_support = pe_support
        self.pe_map_quality = 0

        self.mapq = mapq
        self.insLen = insLen
        self.svt = svt
        self.id = id
        self.homLen = 0

        self.sr_align_quality = 0.0
        self.precise = True
        self.alleles = ""
        self.consensus = ""

    def print_info(self):
        print("Structural Variant Record Information:")
        for key, value in self.__dict__.items():
            print(f"{key}: {value}")


def sort_svs(sv1: StructuralVariantRecord, sv2: StructuralVariantRecord):
    if sv1.chr != sv2.chr:
        return sv1.chr < sv2.chr
    if sv1.sv_start != sv2.sv_start:
        return sv1.sv_start < sv2.sv_start
    if sv1.chr2 != sv2.chr2:
        return sv1.chr2 < sv2.chr2
    if sv1.sv_end != sv2.sv_end:
        return sv1.sv_end < sv2.sv_end
    if sv1.pe_support != sv2.pe_support:
        return sv1.pe_support > sv2.pe_support
    return sv1.sr_support > sv2.sr_support


def is_translocation(svt):
    SVT_TRANS = 5
    return SVT_TRANS <= svt < 9


class EdgeRecord:
    def __init__(self, weight: int, source: int, target: int):
        self.weight = weight
        self.source = source
        self.target = target

    def print_info(self):
        print("EdgeRecord Info:")
        for key, value in self.__dict__.items():
            print(f"{key}: {value}")


TEdgeRecord = EdgeRecord
TEdgeList = List[TEdgeRecord]
TCompEdgeList = Dict[int, TEdgeList]


def sort_edge_records(e1: EdgeRecord, e2: EdgeRecord):
    if e1.weight < e2.weight:
        return -1
    elif e1.weight > e2.weight:
        return 1
    else:
        if e1.source < e2.source:
            return -1
        elif e1.source > e2.source:
            return 1
        else:
            if e1.target < e2.target:
                return -1
            elif e1.target > e2.target:
                return 1
            else:
                return 0


def sv_size_check(start, end, svt):
    # Short reads
    if svt in (0, 1, 2):
        return end - start >= 300
    # Long reads
    elif svt == 3:
        return end - start >= 100
    else:
        return True


def get_sv_type(rec):
    if not rec.flag & 16:  # BAM_FREVERSE
        if not rec.flag & 32:
            return 0  # BAM_FMREVERSE
        else:
            return 2 if rec.pos < rec.mpos else 3
    else:
        if not rec.flag & 32:  # BAM_FMREVERSE
            return 2 if rec.pos > rec.mpos else 3
        else:
            return 1


def hash_string(s):
    h = 37
    for char in s:
        h = (h * 54059) ^ (ord(char) * 76963)
    return h


class Interval:
    def __init__(self, start, end):
        self.start = start
        self.end = end

    @staticmethod
    def right_open(start, end):
        return Interval(start, end - 1)

    def lower(self):
        return self.start

    def upper(self):
        return self.end

    def __hash__(self):
        return hash((self.start, self.end))

    def __eq__(self, other):
        return self.start == other.start and self.end == other.end

    def __contains__(self, item):
        return self.start <= item < self.end

    def __getitem__(self, index):
        if index == 0:
            return self.start
        elif index == 1:
            return self.end
        else:
            raise IndexError("Only indices 0 (for start) and 1 (for end) are valid.")

    def __repr__(self):
        return f"[{self.start}, {self.end})"


class IntervalSet:
    def __init__(self):
        self.intervals = set()

    def insert(self, interval):
        self.intervals.add(interval)

    def __contains__(self, item):
        for interval in self.intervals:
            if item in interval:
                return True
        return False

    def __repr__(self):
        return "{" + ", ".join(map(str, self.intervals)) + "}"

    def __iter__(self):
        return iter(self.intervals)



def getMedian(elements):
    elements = sorted(elements)
    middle = len(elements) // 2
    return elements[middle]


def getPercentile(vec, p):
    idx = int(len(vec) * p)
    return sorted(vec)[idx]