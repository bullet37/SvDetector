import pysam
from utils import *


def bam_name2id(bamfile, ref_name):
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        ref_id = bam.get_tid(ref_name)
        return ref_id


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
                                    print(
                                        "Exclude file needs to be in tab-delimited format (chr, start, end) and start < end.")
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



def merge_sort(pe: List[StructuralVariantRecord], sr: List[StructuralVariantRecord]):
    pe.sort(key=cmp_to_key(sort_svs))
    sr.sort(key=cmp_to_key(sort_svs))
    # Augment PE SVs and append missing SR SVs
    for svt in range(10):
        for i, sr_record in enumerate(sr):
            if sr_record.svt != svt:
                continue
            if sr_record.sr_support == 0 or sr_record.sr_align_quality == 0:
                continue  # SR assembly failed

            # Precise duplicates
            search_window = 500
            sv_exists = False
            for pe_record in pe:
                if abs(pe_record.sv_start - sr_record.sv_start) >= search_window:
                    break
                if pe_record.svt != svt or pe_record.precise:
                    continue
                if pe_record.chr != sr_record.chr or pe_record.chr2 != sr_record.chr2:
                    continue

                # Breakpoints within PE confidence interval?
                if pe_record.sv_start + pe_record.ciposlow < sr_record.sv_start < pe_record.sv_start + pe_record.ciposhigh:
                    if pe_record.sv_end + pe_record.ciendlow < sr_record.sv_end < pe_record.sv_end + pe_record.ciendhigh:
                        sv_exists = True
                        # Augment PE record
                        pe_record.sv_start = sr_record.sv_start
                        pe_record.sv_end = sr_record.sv_end
                        pe_record.ciposlow = sr_record.ciposlow
                        pe_record.ciposhigh = sr_record.ciposhigh
                        pe_record.ciendlow = sr_record.ciendlow
                        pe_record.ciendhigh = sr_record.ciendhigh
                        pe_record.sr_map_quality = sr_record.sr_map_quality
                        pe_record.sr_support = sr_record.sr_support
                        pe_record.insLen = sr_record.insLen
                        pe_record.homLen = sr_record.homLen
                        pe_record.sr_align_quality = sr_record.sr_align_quality
                        pe_record.precise = True
                        pe_record.consensus = sr_record.consensus
                        pe_record.consBp = sr_record.consBp
                        pe_record.mapq += sr_record.mapq

            # SR only SV
            if not sv_exists:
                # Check for PRECISE duplicate
                precise_duplicate = False
                prec_search_window = 10
                for j in range(i + 1, len(sr)):
                    if abs(sr[i].sv_start - sr[j].sv_start) > prec_search_window:
                        break
                    if sr[i].svt != sr[j].svt:
                        continue
                    if sr[i].chr != sr[j].chr or sr[i].chr2 != sr[j].chr2:
                        continue

                    if sr[j].sv_start + sr[j].ciposlow <= sr[i].sv_start <= sr[j].sv_start + sr[j].ciposhigh:
                        if sr[j].sv_end + sr[j].ciendlow <= sr[i].sv_end <= sr[j].sv_end + sr[j].ciendhigh:
                            if sr[i].sr_support < sr[j].sr_support or (i < j and sr[i].sr_support == sr[j].sr_support):
                                precise_duplicate = True

                for j in range(i - 1, -1, -1):
                    if abs(sr[i].sv_start - sr[j].sv_start) > prec_search_window:
                        break
                    if sr[i].svt != sr[j].svt or (sr[i].chr != sr[j].chr) or (sr[i].chr2 != sr[j].chr2):
                        continue  # Mismatching chr

                    # Breakpoints within PE confidence interval?
                    if (sr[j].sv_start + sr[j].ciposlow < sr[i].sv_start) and (
                            sr[i].sv_start < sr[j].sv_start + sr[j].ciposhigh):
                        if (sr[j].sv_end + sr[j].ciendlow < sr[i].sv_end) and (
                                sr[i].sv_end < sr[j].sv_end + sr[j].ciendhigh):
                            # Duplicate, keep better call
                            if (sr[i].sr_support < sr[j].sr_support) or (
                                    (i < j) and (sr[i].sr_support == sr[j].sr_support)):
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
