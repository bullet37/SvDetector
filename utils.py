# Config
from typing import List, Dict, Set
from functools import cmp_to_key
import numpy as np
import pysam
import sys , os
SVT_TRANS = 5

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
        self.outfile = ''
        self.vcffile = ''
        self.genome = ''
        self.exclude = ''
        self.dumpfile = ''
        self.files = []
        self.sampleName = []

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
          if not (os.path.exists(self.exclude) and os.path.isfile(self.exclude) and os.path.getsize(self.exclude) > 0):
              print(f"Exclude file is missing: {self.exclude}", file=sys.stderr)
              sys.exit(1)
          self.hasExcludeFile = True

      if not (os.path.exists(self.genome) and os.path.isfile(self.genome) and os.path.getsize(self.genome) > 0):
        print(f"Reference file is missing: {self.genome}", file=sys.stderr)
        sys.exit(1)

      if self.vcffile:
          if not (os.path.exists(self.vcffile) and os.path.isfile(self.vcffile) and os.path.getsize(self.vcffile) > 0):
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
      if self.minGenoQual < 5:
          self.minGenoQual = 5

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
                      print(f"BAM file chromosome {ref_name} is NOT present in your reference file {self.genome}", file=sys.stderr)
                      sys.exit(1)
          samfile.close()

    def setup(self, args):
      self.genome = args.genome
      self.exclude = args.exclude
      self.outfile = args.outfile
      self.minMapQual = args.map_qual
      self.minTraQual = args.qual_tra
      self.madCutoff = args.mad_cutoff
      self.minClip = args.minclip
      self.minRefSep = args.minrefsep
      self.maxReadSep = args.maxreadsep
      self.vcffile = args.vcffile
      self.minGenoQual = args.geno_qual
      self.dumpfile = args.dump
      self.maxGenoReadCount = args.max_geno_count
      self.files = args.input_files
      self.graphPruning = args.pruning

      self.aliscore = DnaScore(5, -4, -10, -1)
      self.hasDumpFile = bool(args.dump)
      self.minCliqueSize = max(2, args.min_clique_size)
      self.minTraQual = max(self.minMapQual, args.qual_tra)

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
        self.sampleName.append(sample_name)
        samfile.close()

    def _get_SM_tag(self, header, file_name):
        sm_identifiers = set()
        rg_present = False
        lines = header.split('\n')
        for line in lines:
            if line.startswith("@RG"):
                keyvals = line.split("\t")
                for keyval in keyvals:
                    if ":" in keyval:
                        field, value = keyval.split(":", 1)
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

def sort_svs(sv1: StructuralVariantRecord, sv2: StructuralVariantRecord):
    if sv1.chr != sv2.chr:
        return sv1.chr < sv2.chr
    if sv1.svStart != sv2.svStart:
        return sv1.svStart < sv2.svStart
    if sv1.chr2 != sv2.chr2:
        return sv1.chr2 < sv2.chr2
    if sv1.svEnd != sv2.svEnd:
        return sv1.svEnd < sv2.svEnd
    if sv1.peSupport != sv2.peSupport:
        return sv1.peSupport > sv2.peSupport
    return sv1.srSupport > sv2.srSupport


def is_translocation(svt):
    SVT_TRANS = 5
    return SVT_TRANS <= svt < 9

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

def getSVType(rec):
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