import pysam
from utils import *
def bam_name2id(bamfile, ref_name):
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        ref_id = bam.get_tid(ref_name)
        return ref_id

def getLibraryParams(c: Config, validRegions, sampleLib):
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

def parseExcludeIntervals(c: Config, samfile, validRegions):
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

def merge_sort(pe: List[StructuralVariantRecord], sr: List[StructuralVariantRecord]):
    pe.sort(key=cmp_to_key(sort_svs))
    sr.sort(key=cmp_to_key(sort_svs))
    # Augment PE SVs and append missing SR SVs
    for svt in range(10):
        for i, sr_record in enumerate(sr):
            if sr_record.svt != svt:
                continue
            if sr_record.srSupport == 0 or sr_record.srAlignQuality == 0:
                continue  # SR assembly failed

            # Precise duplicates
            search_window = 500
            sv_exists = False
            for pe_record in pe:
                if abs(pe_record.svStart - sr_record.svStart) >= search_window:
                    break
                if pe_record.svt != svt or pe_record.precise:
                    continue
                if pe_record.chr != sr_record.chr or pe_record.chr2 != sr_record.chr2:
                    continue

                # Breakpoints within PE confidence interval?
                if pe_record.svStart + pe_record.ciposlow < sr_record.svStart < pe_record.svStart + pe_record.ciposhigh:
                    if pe_record.svEnd + pe_record.ciendlow < sr_record.svEnd < pe_record.svEnd + pe_record.ciendhigh:
                        sv_exists = True
                        # Augment PE record
                        pe_record.svStart = sr_record.svStart
                        pe_record.svEnd = sr_record.svEnd
                        pe_record.ciposlow = sr_record.ciposlow
                        pe_record.ciposhigh = sr_record.ciposhigh
                        pe_record.ciendlow = sr_record.ciendlow
                        pe_record.ciendhigh = sr_record.ciendhigh
                        pe_record.srMapQuality = sr_record.srMapQuality
                        pe_record.srSupport = sr_record.srSupport
                        pe_record.insLen = sr_record.insLen
                        pe_record.homLen = sr_record.homLen
                        pe_record.srAlignQuality = sr_record.srAlignQuality
                        pe_record.precise = True
                        pe_record.consensus = sr_record.consensus
                        pe_record.consBp = sr_record.consBp
                        pe_record.mapq += sr_record.mapq

            # SR only SV
            if not sv_exists:
                # Check for PRECISE duplicate
                precise_duplicate = False
                prec_search_window = 10
                for j in range(i+1, len(sr)):
                    if abs(sr[i].svStart - sr[j].svStart) > prec_search_window:
                        break
                    if sr[i].svt != sr[j].svt:
                        continue
                    if sr[i].chr != sr[j].chr or sr[i].chr2 != sr[j].chr2:
                        continue

                    if sr[j].svStart + sr[j].ciposlow <= sr[i].svStart <= sr[j].svStart + sr[j].ciposhigh:
                        if sr[j].svEnd + sr[j].ciendlow <= sr[i].svEnd <= sr[j].svEnd + sr[j].ciendhigh:
                            if sr[i].srSupport < sr[j].srSupport or (i < j and sr[i].srSupport == sr[j].srSupport):
                                precise_duplicate = True

                for j in range(i-1, -1, -1):
                    if abs(sr[i].svStart - sr[j].svStart) > prec_search_window:
                        break
                    if sr[i].svt != sr[j].svt or (sr[i].chr != sr[j].chr) or (sr[i].chr2 != sr[j].chr2):
                      continue  # Mismatching chr

                    # Breakpoints within PE confidence interval?
                    if (sr[j].svStart + sr[j].ciposlow < sr[i].svStart) and (sr[i].svStart < sr[j].svStart + sr[j].ciposhigh):
                        if (sr[j].svEnd + sr[j].ciendlow < sr[i].svEnd) and (sr[i].svEnd < sr[j].svEnd + sr[j].ciendhigh):
                            # Duplicate, keep better call
                            if (sr[i].srSupport < sr[j].srSupport) or ((i < j) and (sr[i].srSupport == sr[j].srSupport)):
                                precise_duplicate = True

                if not precise_duplicate:
                    pe.append(sr[i])
                    pe.sort(key=cmp_to_key(sort_svs))



def display_usage():
    print("Usage: SvDetector <command> <arguments>\n")
    print("Short-read SV calling:")
    print("    call         discover and genotype structural variants")
    print("    merge        merge structural variants across VCF/BCF files and within a single VCF/BCF file")
    # print("    filter       filter somatic or germline structural variants\n")
    # print("Long-read SV calling:")
    # print("    lr           long-read SV discovery\n")
    # print("Pan-genome based SV calling (work-in-progress):")
    # print("    pg           pan-genome SV discovery\n")
    # print("Assembly-based SV calling (work-in-progress):")
    # print("    asm          assembly SV site discovery\n")
    # print("Copy-number variant calling:")
    # print("    cnv          discover and genotype copy-number variants")
    # print("    classify     classify somatic or germline copy-number variants\n")
    # print("Multi-sample VCF operations:")
    # print("    markdup      mark duplicate SV sites based on SV allele and GT concordance")
    # print("    compvcf      compare multi-sample VCF file to a ground truth VCF file\n")
    # Deprecated commands are commented out
    # print("Deprecated:")
    # print("    dpe          double paired-end signatures")
    # print("    chimera      ONT chimera flagging\n")