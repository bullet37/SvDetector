import math
import pysam
from datetime import datetime
import numpy as np
from utils import *
from typing import List, Dict, Set


class BoLog:
    def __init__(self):
        smallest_gl = -1000
        self.phred2prob = []
        max_phred = round(-10.0 * smallest_gl)
        for i in range(max_phred + 1):
            self.phred2prob.append(math.pow(10.0, -(i / 10.0)))


def entropy(st):
    alphabet = set(st)
    ent = 0.0
    for c in alphabet:
        ctr = st.count(c)
        freq = ctr / len(st)
        ent += freq * math.log(freq, 2)
    return -ent


def _add_id(svt):
    if svt == 0:
        return "INV"
    elif svt == 1:
        return "INV"
    elif svt == 2:
        return "DEL"
    elif svt == 3:
        return "DUP"
    elif svt == 4:
        return "INS"
    elif svt == 9:
        return "CNV"
    else:
        return "BND"


def _get_span_orientation(svt):
    if is_translocation(svt):
        return svt - 5
    else:
        return svt


def _add_orientation(svt):
    ct = _get_span_orientation(svt)
    if ct == 0:
        return "3to3"
    elif ct == 1:
        return "5to5"
    elif ct == 2:
        return "3to5"
    elif ct == 3:
        return "5to3"
    else:
        return "NtoN"


def _replace_iupac(alleles):
    out = []
    in_tag = 0
    for i in alleles:
        if (in_tag or i in "ACGTNacgtn<>,[]"):
            out.append(i)
            if i == '<':
                in_tag = 1
            elif i == ']':
                in_tag = 2
            elif i == '[':
                in_tag = 3
            elif i == '>' and in_tag == 1:
                in_tag = 0
            elif i == ']' and in_tag == 2:
                in_tag = 0
            elif i == '[' and in_tag == 3:
                in_tag = 0
        else:
            if i in "Uu":
                out.append('T')
            elif i in "Rr":
                out.append('A')
            elif i in "Yy":
                out.append('C')
            elif i in "Ss":
                out.append('C')
            elif i in "Ww":
                out.append('A')
            elif i in "Kk":
                out.append('G')
            elif i in "Mm":
                out.append('A')
            elif i in "Bb":
                out.append('C')
            elif i in "Dd":
                out.append('A')
            elif i in "Hh":
                out.append('A')
            elif i in "Vv":
                out.append('A')
            else:
                out.append('N')
    return ''.join(out)


def _add_alleles(ref, alt):
    return ref + "," + alt


def bcf_gt_is_missing(val):
    return (val >> 1) == 0


def bcf_gt_allele(val):
    return (val >> 1) - 1


def bcf_gt_is_phased(idx):
    return (idx & 1) != 0


def bcf_gt_phased(idx):
    return ((idx + 1) << 1) | 1


def bcf_gt_unphased(idx):
    return ((idx + 1) << 1)


def compute_gls(bl: BoLog, mapq_ref, mapq_alt, gls, gqval, gts, file_c):
    gl = [0.0, 0.0, 0.0]
    bcf_gt_missing = 0
    SMALLEST_GL = -1000
    pe_depth = len(mapq_ref) + len(mapq_alt)
    for mq in mapq_ref:
        prob = bl.phred2prob[mq]
        gl[0] += math.log10(prob)
        gl[1] += math.log10(prob + (1 - prob))
        gl[2] += math.log10(1 - prob)

    for mq in mapq_alt:
        prob = bl.phred2prob[mq]
        gl[0] += math.log10(1 - prob)
        gl[1] += math.log10((1 - prob) + prob)
        gl[2] += math.log10(prob)

    gl[1] -= pe_depth * math.log10(2)

    # Find the best genotype
    gl_best = max(range(3), key=lambda geno: gl[geno])
    gl_best_val = gl[gl_best]
    for geno in range(3):
        gl[geno] = max(gl[geno] - gl_best_val, SMALLEST_GL)

    pl = [int(-10 * gl[geno]) for geno in range(3)]

    if pe_depth and sum(pl) > 0:
        likelihood = math.log10(1 - 1 / sum([bl.phred2prob[p] for p in pl]))
        likelihood = max(likelihood, SMALLEST_GL)
        gqval[file_c] = round(-10 * likelihood)

        if gl_best == 0:
            gts[file_c * 2] = bcf_gt_unphased(1)
            gts[file_c * 2 + 1] = bcf_gt_unphased(1)
        elif gl_best == 1:
            gts[file_c * 2] = bcf_gt_unphased(0)
            gts[file_c * 2 + 1] = bcf_gt_unphased(1)
        else:
            gts[file_c * 2] = bcf_gt_unphased(0)
            gts[file_c * 2 + 1] = bcf_gt_unphased(0)
    else:
        gts[file_c * 2] = bcf_gt_missing
        gts[file_c * 2 + 1] = bcf_gt_missing
        gqval[file_c] = 0

    # update GLs
    gls_index = file_c * 3
    gls[gls_index:gls_index + 3] = [gl[2], gl[1], gl[0]]


# def _add_alleles(ref, chr2, sv, svt):
#     ct = _get_span_orientation(svt)
#     if  is_translocation(svt):
#         if ct == 0:
#             return f"{ref},{ref}]{chr2}:{sv.sv_end}]"
#         elif ct == 1:
#             return f"{ref},[{chr2}:{sv.sv_end}[{ref}"
#         elif ct == 2:
#             return f"{ref},{ref}[{chr2}:{sv.sv_end}["
#         elif ct == 3:
#             return f"{ref},]{chr2}:{sv.sv_end}]{ref}"
#         else:
#             return f"{ref},<{_add_id(svt)}>"
#     else:
#         return f"{ref},<{_add_id(svt)}>"

def vcf_output(c: Config, svs: List[StructuralVariantRecord], jct_count_map, read_count_map, span_count_map):
    samfile = pysam.AlignmentFile(c.files[0], "r")
    bam_hdr = samfile.header
    now = datetime.now()
    bl = BoLog()
    # Output all structural variants
    fmtout = "wb" if c.outfile != "-" else "w"
    fp = pysam.VariantFile(c.outfile, fmtout)
    hdr = fp.header

    # Print VCF header
    hdr.add_meta("fileDate", value=now.strftime("%Y%m%d"))
    hdr.add_meta('ALT', items=[('ID', 'DEL'), ('Description', "Deletion")])
    hdr.add_meta('ALT', items=[('ID', 'DUP'), ('Description', "Duplication")])
    hdr.add_meta('ALT', items=[('ID', 'INV'), ('Description', "Inversion")])
    hdr.add_meta('ALT', items=[('ID', 'BND'), ('Description', "Translocation")])
    hdr.add_meta('ALT', items=[('ID', 'INS'), ('Description', "Insertion")])
    hdr.add_meta('FILTER',
                 items=[('ID', 'LowQual'), ('Description', "Poor quality and insufficient number of PEs and SRs.")])

    hdr.add_meta('INFO', items=[('ID', 'CIEND'), ('Number', 2), ('Type', 'Integer'),
                                ('Description', "PE confidence interval around END")])
    hdr.add_meta('INFO', items=[('ID', 'CIPOS'), ('Number', 2), ('Type', 'Integer'),
                                ('Description', "PE confidence interval around POS")])
    hdr.add_meta('INFO', items=[('ID', 'CHR2'), ('Number', 1), ('Type', 'String'), (
    'Description', "Chromosome for POS2 coordinate in case of an inter-chromosomal translocation")])
    hdr.add_meta('INFO', items=[('ID', 'POS2'), ('Number', 1), ('Type', 'Integer'), (
    'Description', "Genomic position for CHR2 in case of an inter-chromosomal translocation")])
    hdr.add_meta('INFO', items=[('ID', 'END'), ('Number', 1), ('Type', 'Integer'),
                                ('Description', "End position of the structural variant")])
    hdr.add_meta('INFO', items=[('ID', 'PE'), ('Number', 1), ('Type', 'Integer'),
                                ('Description', "Paired-end support of the structural variant")])
    hdr.add_meta('INFO', items=[('ID', 'MAPQ'), ('Number', 1), ('Type', 'Integer'),
                                ('Description', "Median mapping quality of paired-ends")])
    hdr.add_meta('INFO', items=[('ID', 'SRMAPQ'), ('Number', 1), ('Type', 'Integer'),
                                ('Description', "Median mapping quality of split-reads")])
    hdr.add_meta('INFO',
                 items=[('ID', 'SR'), ('Number', 1), ('Type', 'Integer'), ('Description', "Split-read support")])
    hdr.add_meta('INFO', items=[('ID', 'SRQ'), ('Number', 1), ('Type', 'Float'),
                                ('Description', "Split-read consensus alignment quality")])
    hdr.add_meta('INFO', items=[('ID', 'CONSENSUS'), ('Number', 1), ('Type', 'String'),
                                ('Description', "Split-read consensus sequence")])
    hdr.add_meta('INFO', items=[('ID', 'CONSBP'), ('Number', 1), ('Type', 'Integer'),
                                ('Description', "Consensus SV breakpoint position")])
    hdr.add_meta('INFO',
                 items=[('ID', 'CE'), ('Number', 1), ('Type', 'Float'), ('Description', "Consensus sequence entropy")])
    hdr.add_meta('INFO', items=[('ID', 'CT'), ('Number', 1), ('Type', 'String'),
                                ('Description', "Paired-end signature induced connection type")])
    hdr.add_meta('INFO', items=[('ID', 'SVLEN'), ('Number', 1), ('Type', 'Integer'),
                                ('Description', "Insertion length for SVTYPE=INS.")])
    hdr.add_meta('INFO', items=[('ID', 'IMPRECISE'), ('Number', 0), ('Type', 'Flag'),
                                ('Description', "Imprecise structural variation")])
    hdr.add_meta('INFO', items=[('ID', 'PRECISE'), ('Number', 0), ('Type', 'Flag'),
                                ('Description', "Precise structural variation")])
    hdr.add_meta('INFO', items=[('ID', 'SVTYPE'), ('Number', 1), ('Type', 'String'),
                                ('Description', "Type of structural variant")])
    hdr.add_meta('INFO', items=[('ID', 'SVMETHOD'), ('Number', 1), ('Type', 'String'),
                                ('Description', "Type of approach used to detect SV")])
    hdr.add_meta('INFO', items=[('ID', 'INSLEN'), ('Number', 1), ('Type', 'Integer'),
                                ('Description', "Predicted length of the insertion")])
    hdr.add_meta('INFO', items=[('ID', 'HOMLEN'), ('Number', 1), ('Type', 'Integer'),
                                ('Description', "Predicted microhomology length using a max. edit distance of 2")])
    hdr.add_meta('FORMAT', items=[('ID', 'GT'), ('Number', 1), ('Type', 'String'), ('Description', "Genotype")])
    hdr.add_meta('FORMAT', items=[('ID', 'GL'), ('Number', 3), ('Type', 'Float'),
                                  ('Description', "Log10-scaled genotype likelihoods for RR,RA,AA genotypes")])
    hdr.add_meta('FORMAT',
                 items=[('ID', 'GQ'), ('Number', 1), ('Type', 'Integer'), ('Description', "Genotype Quality")])
    hdr.add_meta('FORMAT',
                 items=[('ID', 'FT'), ('Number', 1), ('Type', 'String'), ('Description', "Per-sample genotype filter")])
    hdr.add_meta('FORMAT', items=[('ID', 'RC'), ('Number', 1), ('Type', 'Integer'),
                                  ('Description', "Raw high-quality read counts or base counts for the SV")])
    hdr.add_meta('FORMAT', items=[('ID', 'RCL'), ('Number', 1), ('Type', 'Integer'), (
    'Description', "Raw high-quality read counts or base counts for the left control region")])
    hdr.add_meta('FORMAT', items=[('ID', 'RCR'), ('Number', 1), ('Type', 'Integer'), (
    'Description', "Raw high-quality read counts or base counts for the right control region")])
    hdr.add_meta('FORMAT', items=[('ID', 'RDCN'), ('Number', 1), ('Type', 'Integer'),
                                  ('Description', "Read-depth based copy-number estimate for autosomal sites")])
    hdr.add_meta('FORMAT', items=[('ID', 'DR'), ('Number', 1), ('Type', 'Integer'),
                                  ('Description', "# high-quality reference pairs")])
    hdr.add_meta('FORMAT', items=[('ID', 'DV'), ('Number', 1), ('Type', 'Integer'),
                                  ('Description', "# high-quality variant pairs")])
    hdr.add_meta('FORMAT', items=[('ID', 'RR'), ('Number', 1), ('Type', 'Integer'),
                                  ('Description', "# high-quality reference junction reads")])
    hdr.add_meta('FORMAT', items=[('ID', 'RV'), ('Number', 1), ('Type', 'Integer'),
                                  ('Description', "# high-quality variant junction reads")])

    # Add reference
    hdr.add_meta("reference", value=str(c.genome))
    for i in range(len(bam_hdr.references)):
        hdr.contigs.add(bam_hdr.references[i], length=bam_hdr.lengths[i])

    # Add samples
    for file_c in range(len(c.files)):
        hdr.add_sample(c.sampleName[file_c])

    if svs:
        num_samples = len(hdr.samples)
        gts = np.zeros(num_samples * 2, dtype=np.int32)
        gls = np.zeros(num_samples * 3, dtype=np.float32)
        rcl = np.zeros(num_samples, dtype=np.int32)
        rc = np.zeros(num_samples, dtype=np.int32)
        rcr = np.zeros(num_samples, dtype=np.int32)
        cnest = np.zeros(num_samples, dtype=np.int32)
        drcount = np.zeros(num_samples, dtype=np.int32)
        dvcount = np.zeros(num_samples, dtype=np.int32)
        rrcount = np.zeros(num_samples, dtype=np.int32)
        rvcount = np.zeros(num_samples, dtype=np.int32)
        gqval = np.zeros(num_samples, dtype=np.int32)
        ftarr = [''] * num_samples

        # Iterate all structural variants
        for sv in svs:
            if sv.sr_support == 0 and sv.pe_support == 0:
                continue

            # In discovery mode, skip SVs that have less than 2 reads support after genotyping
            if not c.hasVcfFile:
                total_gt_sup = sum(
                    [len(span_count_map[file_c][sv.id].alt) + len(jct_count_map[file_c][sv.id].alt) for file_c in
                     range(len(c.files))])
                if total_gt_sup < 2:
                    continue

            # Output main vcf fields, 输出 VCF 主要字段
            rec = hdr.new_record()
            rec.filter.add('PASS')
            if sv.chr == sv.chr2:
                if sv.pe_support < 3 or sv.pe_map_quality < 20 or sv.sr_support < 3 or sv.sr_map_quality < 20:
                    rec.filter.add('LowQual')
            else:
                if sv.pe_support < 5 or sv.pe_map_quality < 20 or sv.sr_support < 5 or sv.sr_map_quality < 20:
                    rec.filter.add('LowQual')

            rec.chrom = bam_hdr.get_reference_name(sv.chr)
            sv_start_pos = max(sv.sv_start - 1, 1)

            chrom_name = bam_hdr.references[sv.chr2]
            sv_end_pos = sv.sv_end if sv.chr == sv.chr2 else sv.sv_end - 1
            sv_end_pos = max(min(sv_end_pos, bam_hdr.get_reference_length(chrom_name) - 1), 1)
            rec.pos = sv_start_pos
            id = _add_id(sv.svt) + str(sv.id).zfill(8)
            rec.id = id

            rec.alleles = ["A", "T"]
            # todo :  _generateProbes
            # alleles = _replace_iupac(sv.alleles)
            # rec.alleles = [alleles]

            # Add INFO fields
            if sv.precise:
                rec.info['PRECISE'] = True
            else:
                rec.info['IMPRECISE'] = True

            rec.info['SVTYPE'] = _add_id(sv.svt)
            SvDetector_version_number = "1.1.8"
            SvDetector_version = "EMBL.SvDetector" + SvDetector_version_number
            rec.info['SVMETHOD'] = SvDetector_version

            if sv.svt < 5:
                rec.stop = sv_end_pos
            else:
                rec.stop = sv_start_pos + 2
                rec.info['CHR2'] = hdr.get_reference_name(sv.chr2)
                rec.info['POS2'] = sv_end_pos

            if sv.svt == 4:
                rec.info['SVLEN'] = sv.insLen

            rec.info['PE'] = sv.pe_support
            rec.info['MAPQ'] = sv.pe_map_quality
            rec.info['CT'] = _add_orientation(sv.svt)
            ciend = [sv.ciendlow, sv.ciendhigh]
            cipos = [sv.ciposlow, sv.ciposhigh]
            rec.info['CIPOS'] = cipos
            rec.info['CIEND'] = ciend

            if sv.precise:
                rec.info['SRMAPQ'] = sv.sr_map_quality
                rec.info['INSLEN'] = sv.insLen
                rec.info['HOMLEN'] = sv.homLen
                rec.info['SR'] = sv.sr_support
                rec.info['SRQ'] = sv.sr_align_quality
                if sv.consensus:
                    rec.info['CONSENSUS'] = sv.consensus
                    rec.info['CE'] = entropy(sv.consensus)
                    rec.info['CONSBP'] = sv.cons_bp

            # Add genotype columns
            for file_c in range(len(c.files)):
                rcl[file_c] = 0
                rc[file_c] = 0
                rcr[file_c] = 0
                cnest[file_c] = 0
                drcount[file_c] = 0
                dvcount[file_c] = 0
                rrcount[file_c] = 0
                rvcount[file_c] = 0
                drcount[file_c] = len(span_count_map[file_c][sv.id].ref)
                dvcount[file_c] = len(span_count_map[file_c][sv.id].alt)
                rrcount[file_c] = len(jct_count_map[file_c][sv.id].ref)
                rvcount[file_c] = len(jct_count_map[file_c][sv.id].alt)

                #  Compute GLs
                if sv.precise:
                    compute_gls(bl, jct_count_map[file_c][sv.id].ref, jct_count_map[file_c][sv.id].alt, gls, gqval, gts,
                                file_c)
                else:
                    compute_gls(bl, span_count_map[file_c][sv.id].ref, span_count_map[file_c][sv.id].alt, gls, gqval,
                                gts, file_c)

                # 更新 RCs
                rcl[file_c] = read_count_map[file_c][sv.id].leftRC
                rc[file_c] = read_count_map[file_c][sv.id].rc
                rcr[file_c] = read_count_map[file_c][sv.id].rightRC
                cnest[file_c] = -1
                if rcl[file_c] + rcr[file_c] > 0:
                    cnest[file_c] = round(2 * rc[file_c] / (rcl[file_c] + rcr[file_c]))

                # Genotype filter
                ftarr[file_c] = "LowQual" if gqval[file_c] < 15 else "PASS"

            qvalout = max(0, min(sv.mapq, 10000))
            rec.qual = qvalout

            for i, sample in enumerate(hdr.samples):
                gl_values = tuple([float(gls[i * 3 + j]) for j in range(3)])
                rec.samples[sample]['GL'] = gl_values
                rec.samples[sample]['GQ'] = int(gqval[i])
                rec.samples[sample]['FT'] = ftarr[i]
                rec.samples[sample]['RCL'] = int(rcl[i])
                rec.samples[sample]['RC'] = int(rc[i])
                rec.samples[sample]['RCR'] = int(rcr[i])
                rec.samples[sample]['RDCN'] = int(cnest[i])
                rec.samples[sample]['DR'] = int(drcount[i])
                rec.samples[sample]['DV'] = int(dvcount[i])
                rec.samples[sample]['RR'] = int(rrcount[i])
                rec.samples[sample]['RV'] = int(rvcount[i])

            fp.write(rec)
    samfile.close()
    fp.close()
    if c.outfile != "-":
        pysam.tabix_index(c.outfile, preset="vcf", force=True)
