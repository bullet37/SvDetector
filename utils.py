import pathlib
from typing import List, Dict, Set
from functools import cmp_to_key
import hashlib
import heapq
import copy
import datetime
import numpy as np

DELLY_SVT_TRANS = 5

class Config:
    def __init__(self):
        self.minMapQual = 1
        self.minTraQual = 20
        self.minGenoQual = 5
        self.madCutoff = 9
        self.madNormalCutoff = 5
        self.nchr = 0
        self.minimumFlankSize = 13
        self.indelsize = 1000
        self.minConsWindow = 100
        self.graphPruning = 1000
        self.minRefSep = 25
        self.maxReadSep = 40
        self.minClip = 25
        self.maxGenoReadCount = 250
        self.minCliqueSize = 2
        self.flankQuality = 0.95
        self.hasExcludeFile = False
        self.hasVcfFile = False
        self.hasDumpFile = False
        self.svtset = set()
        self.aliscore = DnaScore(5, -4, -10, -1)
        self.outfile = pathlib.Path('')
        self.vcffile = pathlib.Path('')
        self.genome = pathlib.Path('')
        self.exclude = pathlib.Path('')
        self.dumpfile = pathlib.Path('')
        self.files = []
        self.sampleName = []

    def print_info(self):
        print("Config Information:")
        for key, value in self.__dict__.items():
            print(f"{key}: {value}")


class DnaScore:
    def __init__(self, m=5, mm=-4, gapopen=-10, gapextension=-1):
        self.match = m
        self.mismatch = mm
        self.go = gapopen
        self.ge = gapextension
        self.inf = 1000000

    def print_info(self):
        print("DnaScore Information:")
        for key, value in self.__dict__.items():
            print(f"{key}: {value}")


class StructuralVariantRecord:
    def __init__(self, chr=0, sv_start=0, chr2=0, svEnd=0, ciposlow=0, ciposhigh=0, ciendlow=0, ciendhigh=0, srSupport=0, peSupport=0,srMapQuality=0, mapq=0, insLen=0, svt=-1, id=0):
        self.chr = chr
        self.svStart = sv_start
        self.chr2 = chr2
        self.svEnd = svEnd
        self.ciposlow = ciposlow
        self.ciposhigh = ciposhigh
        self.ciendlow = ciendlow
        self.ciendhigh = ciendhigh
        self.srSupport = srSupport
        self.srMapQuality = srMapQuality
        self.mapq = mapq
        self.insLen = insLen
        self.svt = svt
        self.id = id
        self.homLen = 0
        self.peSupport = peSupport
        self.peMapQuality = 0
        self.srAlignQuality = 0.0
        self.precise = True
        self.alleles = ""
        self.consensus = ""

    def print_info(self):
        print("Structural Variant Record Information:")
        for key, value in self.__dict__.items():
            print(f"{key}: {value}")

class LibraryInfo:
    def __init__(self, rs=0, median=0, mad=0, minNormalISize=0,
                 minISizeCutoff=0, maxNormalISize=0, maxISizeCutoff=0, abnormal_pairs=0):
        self.rs = rs
        self.median = 0
        self.mad = 0
        self.minNormalISize = 0
        self.minISizeCutoff = 0
        self.maxNormalISize = 0
        self.maxISizeCutoff = 0
        self.abnormal_pairs = 0

    def print_info(self):
        print("Library Info:")
        for key, value in self.__dict__.items():
            print(f"{key}: {value}")



def is_translocation(svt):
    DELLY_SVT_TRANS = 5
    return DELLY_SVT_TRANS <= svt < 9

class EdgeRecord:
    def __init__(self,  weight: int,  source: int,  target: int):
        self.weight =  weight
        self.source =  source
        self.target =  target

    def print_info(self):
        print("EdgeRecord Info:")
        for key, value in self.__dict__.items():
            print(f"{key}: {value}")

TEdgeRecord = EdgeRecord
TEdgeList = List[TEdgeRecord]
TCompEdgeList = Dict[int, TEdgeList]

def sort_edge_records(e1, e2):
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

class Interval:
    def __init__(self, start, end):
        self.start = start
        self.end = end
    @staticmethod
    def right_open(start, end):
        return Interval(start, end-1)

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


import pysam

def bam_name2id(bamfile, ref_name):
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        ref_id = bam.get_tid(ref_name)
        return ref_id

def getLibraryParams(c, validRegions, sampleLib):
    samfiles = [pysam.AlignmentFile(f, 'r') for f in c.files]
    for f in c.files:
        pysam.index(f)
    for file_c, file_name in enumerate(c.files):
        maxAlignmentsScreened = 10000000
        maxNumAlignments = 1000000
        minNumAlignments = 1000
        alignmentCount = 0
        processedNumPairs = 0
        processedNumReads = 0
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
                        if alignmentCount > maxAlignmentsScreened or (processedNumReads >= maxNumAlignments and processedNumPairs == 0) or processedNumPairs >= maxNumAlignments:
                            libCharacterized = True
                            break
                        alignmentCount += 1

                        if processedNumReads < maxNumAlignments:
                            readSize.append(rec.infer_read_length())
                            processedNumReads += 1

                        if rec.is_paired and not rec.mate_is_unmapped and rec.rname == rec.mrnm:
                            if processedNumPairs < maxNumAlignments:
                                vecISize.append(abs(rec.tlen))
                                if getSVType(rec) == 2:
                                    rplus += 1
                                else:
                                    nonrplus += 1
                                processedNumPairs += 1
                if libCharacterized:
                    break
            if libCharacterized:
                break

        if processedNumReads >= minNumAlignments:
            readSize.sort()
            sampleLib[file_c].rs = readSize[len(readSize) // 2] # Median

        if processedNumPairs >= minNumAlignments:
            vecISize.sort()
            median = vecISize[len(vecISize) // 2]
            absDev = [abs(i - median) for i in vecISize]
            absDev.sort()
            mad = absDev[len(absDev) // 2]

            if 50 <= median <= 100000:
                if rplus < nonrplus:
                    print(f"Warning: Sample has a non-default paired-end layout! File: {file_name}")
                    print("The expected paired-end orientation is ---Read1---> <---Read2--- which is the default illumina paired-end layout.")
                else:
                    lib = sampleLib[file_c]
                    lib.median = median
                    lib.mad = mad  # Median Absolute Deviation
                    lib.maxNormalISize = median + (c.madNormalCutoff * mad)  # 正常范围threshold
                    lib.minNormalISize = max(median - (c.madNormalCutoff * mad), 0)
                    lib.maxISizeCutoff = max(median + (c.madCutoff * mad), max(2 * lib.rs, 500)) # 一个异常插入片段的阈值
                    lib.minISizeCutoff = max(median - (c.madCutoff * mad), 0)

def _parseExcludeIntervals(c, samfile, validRegions):
    n_targets = len(samfile.references)
    exclg = [set() for _ in range(n_targets)]
    valid_chr = [True] * n_targets

    if c.hasExcludeFile:
        with open(c.exclude, 'r') as chr_file:
            for line in chr_file:
                tokens = line.split()
                if tokens:
                    chr_name = tokens[0]
                    tid = bam_name2id(samfile, chr_name)
                    if tid >= 0:
                        if len(tokens) > 1:
                            try:
                                start = int(tokens[1])
                            except ValueError:
                                print("Exclude file needs to be in tab-delimited format: chr, start, end")
                                print(f"Offending line: {line}")
                                return False
                            if len(tokens) > 2:
                                end = int(tokens[2])
                                if start < end:
                                    exclg[tid].add(Interval.right_open(start, end))
                                else:
                                    print("Exclude file needs to be in tab-delimited format (chr, start, end) and start < end.")
                                    print(f"Offending line: {line}")
                                    return False
                            else:
                                print("Exclude file needs to be in tab-delimited format: chr, start, end")
                                print(f"Offending line: {line}")
                                return False
                        else:
                            valid_chr[tid] = False

    for i in range(n_targets):
        if not valid_chr[i]:
            continue
        istart = 0
        for it in sorted(exclg[i], key=lambda x: x.start):
            if istart + 1 < it.lower():
                validRegions[i].insert(Interval.right_open(istart, it.lower() - 1))
            istart = it.upper()
        if istart + 1 < samfile.lengths[i]:
            validRegions[i].insert(Interval.right_open(istart, samfile.lengths[i]))
    return True


def getMedian(elements):
    elements = sorted(elements)
    middle = len(elements) // 2
    return elements[middle]


def getPercentile(vec, p):
    idx = int(len(vec) * p)
    return sorted(vec)[idx]