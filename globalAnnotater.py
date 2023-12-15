# # todo
# class BpRegion:
#     def __init__(self, regionStart=0, regionEnd=0, bppos=0, homLeft=0, homRight=0, svt=0, id=0, bpPoint=0):
#         self.regionStart = regionStart
#         self.regionEnd = regionEnd
#         self.bppos = bppos
#         self.homLeft = homLeft
#         self.homRight = homRight
#         self.svt = svt
#         self.id = id
#         self.bpPoint = bpPoint
#
# class SpanPoint:
#     def __init__(self, bppos=0, svt=0, id=0, chr2=0, otherBppos=0):
#         self.bppos = bppos
#         self.svt = svt
#         self.id = id
#         self.chr2 = chr2
#         self.otherBppos = otherBppos
#
# class FileCoverage:
#     def __init__(self, target_len):
#         self.cov_fragment = [0] * target_len
#         self.cov_bases = [0] * target_len
#         self.bp_occupied = bitarray(target_len)
#         self.bp_occupied.setall(False)
#         self.span_bp = bitarray(target_len)
#         self.span_bp.setall(False)
#         self.span_points = []
#
# class JunctionCount:
#     def __init__(self):
#         self.ref = []
#         self.alt = []
#
# class ReadCount:
#     def __init__(self, l=0, m=0, r=0):
#         self.leftRC = l
#         self.rc = m
#         self.rightRC = r
#
# class SpanningCount:
#     def __init__(self):
#         self.ref = []
#         self.alt = []
#
#
# import datetime
# import gzip
# import pysam
# from collections import defaultdict
# from itertools import groupby
# from operator import itemgetter
# from bitarray import bitarray
#
#
# def _generate_probes(config, header, svs, refProbeArr, consProbeArr, bpRegion, svOnChr):
#     # Implement probe generation logic here...
#     pass
#
# # annotateCoverage(c, sampleLib, svs, rcMap, jctMap, spanMap);
# def annotate_coverage(config, sampleLib, svs, covCount, countMap, spanMap):
#     # Open file handles
#     samfiles = []
#     headers = []
#     totalTarget = 0
#     for file_path in config.files:
#         samfile = pysam.AlignmentFile(file_path, "r")
#         samfiles.append(samfile)
#         header = samfile.header
#         headers.append(header)
#         totalTarget += len(header.references)
#
#     # Initialize coverage count maps
#     for file_c in range(len(config.files)):
#         covCount.append([0] * len(svs))
#         countMap.append([0] * len(svs))
#         spanMap.append([0] * len(svs))
#
#     # Reference and consensus probes
#     refProbeArr = [[None] * len(svs), [None] * len(svs)]
#     consProbeArr = [[None] * len(svs), [None] * len(svs)]
#     bpRegion = [[] for _ in range(len(headers[0].references))]
#     sv_on_chr = [False] * len(headers[0].references)
#
#     # Generate probes
#     _generate_probes(config, headers[0], svs, refProbeArr, consProbeArr, bpRegion, sv_on_chr)
#
#     # Initialize counts
#     refAlignedReadCount = [[0] * len(svs) for _ in range(len(config.files))]
#     refAlignedSpanCount = [[0] * len(svs) for _ in range(len(config.files))]
#
#     # Dump file
#     if config.hasDumpFile:
#         with gzip.open(config.dumpfile, 'wt') as dumpOut:
#             dumpOut.write("#svid\tbam\tqname\tchr\tpos\tmatechr\tmatepos\tmapq\ttype\n")
#
#
#     file_coverage = []
#
#     for file_c, sam_file in enumerate(samfiles):
#         qualities = {}
#         qualitiestra = {}
#         clip = {}
#         cliptra = {}
#
#         # Iterate chromosomes
#         for ref_index in range(len(headers[file_c].references)):
#             if not sv_on_chr(ref_index, svs):
#                 continue
#
#             nodata = True
#             file_path = config.files[file_c]
#             if file_path.endswith("cram"):
#                 nodata = False
#
#             mapped = sam_file.get_index_statistics()[ref_index].mapped
#             if mapped:
#                 nodata = False
#             if nodata:
#                 continue
#
#             # Coverage track
#             max_coverage = 2**16 - 1
#             cov_fragment = [0] * hdr[file_c].target_len[refIndex]
#             cov_bases = [0] * hdr[file_c].target_len[refIndex]
#
#             # Flag breakpoint regions
#             bp_occupied = [0] * hdr[file_c].target_len[refIndex]
#             for i in range(len(bpRegion[refIndex])):
#                 for k in range(bpRegion[refIndex][i].regionStart, bpRegion[refIndex][i].regionEnd):
#                     bp_occupied[k] = 1