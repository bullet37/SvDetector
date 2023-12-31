import pysam
from utils import *


def is_translocation(svt):
    SVT_TRANS = 5
    return SVT_TRANS <= svt < 9


def is_rec_translocation(rec):
    return rec.reference_id != rec.next_reference_id


class PEAlignRecord:
    def __init__(self, rec, pair_quality, alen, malen, median, mad, max_i_size):
        self.tid = rec.reference_id
        self.pos = rec.reference_start
        self.mtid = rec.next_reference_id
        self.mpos = rec.next_reference_start
        self.alen = alen
        self.malen = malen
        self.median = median
        self.mad = mad
        self.max_normal_isize = max_i_size
        self.flag = rec.flag
        self.map_quality = pair_quality

    def print_info(self):
        print("PEAlignRecord Info")
        for key, value in self.__dict__.items():
            print(f"{key}: {value}")


def sort_PERecords(s1, s2):
    if s1.tid == s1.mtid:
        if min(s1.pos, s1.mpos) != min(s2.pos, s2.mpos):
            return min(s1.pos, s1.mpos) - min(s2.pos, s2.mpos)
        if max(s1.pos, s1.mpos) != max(s2.pos, s2.mpos):
            return max(s1.pos, s1.mpos) - max(s2.pos, s2.mpos)
        return s1.max_normal_isize - s2.max_normal_isize
    else:
        if s1.pos != s2.pos:
            return s1.pos - s2.pos
        if s1.mpos != s2.mpos:
            return s1.mpos - s2.mpos
        return s1.max_normal_isize - s2.max_normal_isize


def get_span_orientation(svt: int):
    if is_translocation(svt):
        return svt - 5
    else:
        return svt


def init_clique(el: PEAlignRecord, sv_start: int, sv_end: int, wiggle: int, svt: int):
    if is_translocation(svt):
        ct = get_span_orientation(svt)
        if ct % 2 == 0:
            sv_start = el.pos + el.alen
            sv_end = el.mpos if ct >= 2 else el.mpos + el.malen
        else:
            sv_start = el.pos
            sv_end = el.mpos + el.malen if ct >= 2 else el.mpos
        wiggle = el.max_normal_isize
    else:
        if svt == 0:
            sv_start = el.mpos + el.malen
            sv_end = el.pos + el.alen
            wiggle = el.max_normal_isize - max(el.alen, el.malen)
        elif svt == 1:
            sv_start = el.mpos
            sv_end = el.pos
            wiggle = el.max_normal_isize - max(el.alen, el.malen)
        elif svt == 2:
            sv_start = el.mpos + el.malen
            sv_end = el.pos
            wiggle = -el.max_normal_isize
        elif svt == 3:
            sv_start = el.mpos
            sv_end = el.pos + el.alen
            wiggle = el.max_normal_isize
    return sv_start, sv_end, wiggle


def update_clique(el: PEAlignRecord, sv_start: int, sv_end: int, wiggle: int, svt: int):
    if is_translocation(svt):
        ct = get_span_orientation(svt)
        new_sv_start = new_sv_end = new_wiggle = None
        if ct % 2 == 0:
            new_sv_start = max(sv_start, el.pos + el.alen)
            new_wiggle = wiggle - (new_sv_start - sv_start)
            if ct >= 2:
                new_sv_end = min(sv_end, el.mpos)
            else:
                new_sv_end = max(sv_end, el.mpos + el.malen)
            new_wiggle -= (sv_end - new_sv_end)
        else:
            new_sv_start = min(sv_start, el.pos)
            new_wiggle = wiggle - (sv_start - new_sv_start)
            if ct >= 2:
                new_sv_end = max(sv_end, el.mpos + el.malen)
            else:
                new_sv_end = min(sv_end, el.mpos)
            new_wiggle -= (sv_end - new_sv_end)

        if new_wiggle > 0:
            sv_start = new_sv_start
            sv_end = new_sv_end
            wiggle = new_wiggle
            return True, sv_start, sv_end, wiggle
        return False, sv_start, sv_end, wiggle
    else:
        if svt in [0, 1]:
            ct = get_span_orientation(svt)
            new_sv_start = new_sv_end = new_wiggle = None
            wiggle_change = None
            if ct == 0:
                new_sv_start = max(sv_start, el.mpos + el.malen)
                new_sv_end = max(sv_end, el.pos + el.alen)
                new_wiggle = min(el.max_normal_isize - (new_sv_start - el.mpos),
                                 el.max_normal_isize - (new_sv_end - el.pos))
                wiggle_change = wiggle - max(new_sv_start - sv_start, new_sv_end - sv_end)
            else:
                new_sv_start = min(sv_start, el.mpos)
                new_sv_end = min(sv_end, el.pos)
                new_wiggle = min(el.max_normal_isize - (el.mpos + el.malen - new_sv_start),
                                 el.max_normal_isize - (el.pos + el.alen - new_sv_end))
                wiggle_change = wiggle - max(sv_start - new_sv_start, sv_end - new_sv_end)

            if wiggle_change < new_wiggle:
                new_wiggle = wiggle_change

            if new_sv_start < new_sv_end and new_wiggle >= 0:
                sv_start = new_sv_start
                sv_end = new_sv_end
                wiggle = new_wiggle
                return True, sv_start, sv_end, wiggle

        elif svt == 2:
            new_sv_start = max(sv_start, el.mpos + el.malen)
            new_sv_end = min(sv_end, el.pos)
            new_wiggle = el.pos + el.alen - el.mpos - el.max_normal_isize - (new_sv_end - new_sv_start)
            wiggle_change = wiggle + (sv_end - sv_start) - (new_sv_end - new_sv_start)

            if wiggle_change > new_wiggle:
                new_wiggle = wiggle_change

            if new_sv_start < new_sv_end and new_wiggle <= 0:
                sv_start = new_sv_start
                sv_end = new_sv_end
                wiggle = new_wiggle
                return True, sv_start, sv_end, wiggle

        elif svt == 3:
            new_sv_start = min(sv_start, el.mpos)
            new_sv_end = max(sv_end, el.pos + el.alen)
            new_wiggle = el.pos - (el.mpos + el.malen) + el.max_normal_isize - (new_sv_end - new_sv_start)
            wiggle_change = wiggle - ((new_sv_end - new_sv_start) - (sv_end - sv_start))

            if wiggle_change < new_wiggle:
                new_wiggle = wiggle_change

            if new_sv_start < new_sv_end and new_wiggle >= 0:
                sv_start = new_sv_start
                sv_end = new_sv_end
                wiggle = new_wiggle
                return True, sv_start, sv_end, wiggle
    return False, sv_start, sv_end, wiggle


def min_coord_fn(position: int, mpos: int, svt: int):
    if is_translocation(svt): return position
    return min(position, mpos)


def max_coord_fn(position: int, mpos: int, svt: int):
    if is_translocation(svt): return mpos
    return max(position, mpos)


def hash_combine(seed, value):
    return seed ^ (hash(value) + 0x9e3779b9 + (seed << 6) + (seed >> 2))


def hash_pair(rec):
    seed = hash_string(rec.qname)
    seed = hash_combine(seed, rec.tid)
    seed = hash_combine(seed, rec.pos)
    seed = hash_combine(seed, rec.next_reference_id)
    seed = hash_combine(seed, rec.mpos)
    return seed


def hash_pair_mate(rec):
    seed = hash_string(rec.qname)
    seed = hash_combine(seed, rec.next_reference_id)
    seed = hash_combine(seed, rec.mpos)
    seed = hash_combine(seed, rec.tid)
    seed = hash_combine(seed, rec.pos)
    return seed


def first_pair_obs(rec, last_aligned_pos_reads):
    if rec.reference_id == rec.next_reference_id:
        return (rec.reference_start < rec.next_reference_start) or \
            ((rec.reference_start == rec.next_reference_start) and (
                        hash_string(rec.query_name) not in last_aligned_pos_reads))
    else:
        return rec.reference_id < rec.next_reference_id


def alignmentLength(rec):
    cigar = rec.cigartuples
    alen = 0
    for op, length in cigar:
        if op in [pysam.CMATCH, pysam.CEQUAL, pysam.CDIFF, pysam.CDEL, pysam.CREF_SKIP]:
            alen += length
    return alen


def sv_size_check(s, e, svt, inslen=None):
    diff = e - s
    if svt == 0 or svt == 1 or svt == 2:
        if inslen:
            return diff >= 15
        return diff >= 300
    elif svt == 3:
        if inslen:
            return diff >= 15
        return diff >= 100
    elif svt == 4:
        return inslen >= 15 if inslen else True
    return True


def isizeMappingPos(rec, max_isize_cutoff):
    if is_rec_translocation(rec):
        orient = get_sv_type(rec)
        return 4 + orient
    else:
        if rec.pos == rec.mpos:
            return -1  # No SV
        orient = get_sv_type(rec)
        if orient == 0:
            return 0
        elif orient == 1:
            return 1
        elif orient == 2:
            if max_isize_cutoff > abs(rec.isize):
                return -1
            else:
                return 2
        else:
            if abs(rec.pos - rec.mpos) < 100:
                return -1  # Too small
            return 3


def pairs_disagree(min_coord, max_coord, alen, max_normal_isize,
                   min_coord_next, max_coord_next, alen_next, max_normal_isize_next, svt):
    if is_translocation(svt):
        ct = get_span_orientation(svt)
        if ct % 2 == 0:
            if (min_coord_next + alen_next - min_coord) > max_normal_isize:
                return True
            if ct >= 2:
                if max_coord_next < max_coord:
                    if (max_coord + alen - max_coord_next) > max_normal_isize:
                        return True
                else:
                    if (max_coord_next + alen_next - max_coord) > max_normal_isize_next:
                        return True
            else:
                if max_coord_next < max_coord:
                    if (max_coord + alen - max_coord_next) > max_normal_isize_next:
                        return True
                else:
                    if (max_coord_next + alen_next - max_coord) > max_normal_isize:
                        return True
        else:
            if (min_coord_next + alen_next - min_coord) > max_normal_isize_next:
                return True
            if ct >= 2:
                if max_coord_next < max_coord:
                    if (max_coord + alen - max_coord_next) > max_normal_isize_next:
                        return True
                else:
                    if (max_coord_next + alen_next - max_coord) > max_normal_isize:
                        return True
            else:
                if max_coord_next < max_coord:
                    if (max_coord + alen - max_coord_next) > max_normal_isize_next:
                        return True
                else:
                    if (max_coord_next + alen_next - max_coord) > max_normal_isize:
                        return True
        return False
    else:
        if svt < 2:
            if not svt:
                # 左跨越倒置
                if (min_coord_next + alen_next - min_coord) > max_normal_isize:
                    return True
                if max_coord_next < max_coord and (max_coord + alen - max_coord_next) > max_normal_isize_next:
                    return True
                if max_coord_next >= max_coord and (max_coord_next + alen_next - max_coord) > max_normal_isize:
                    return True
            else:
                # 右跨越倒置
                if (min_coord_next + alen_next - min_coord) > max_normal_isize_next:
                    return True
                if max_coord_next < max_coord and (max_coord + alen - max_coord_next) > max_normal_isize:
                    return True
                if max_coord_next >= max_coord and (max_coord_next + alen_next - max_coord) > max_normal_isize_next:
                    return True
            return False
        elif svt == 2:
            # 缺失
            if (min_coord_next + alen_next - min_coord) > max_normal_isize:
                return True
            if max_coord_next < max_coord and (max_coord + alen - max_coord_next) > max_normal_isize:
                return True
            if max_coord_next >= max_coord and (max_coord_next + alen_next - max_coord) > max_normal_isize_next:
                return True
            if max_coord < min_coord_next or max_coord_next < min_coord:
                return True
            return False
        elif svt == 3:
            if (min_coord_next + alen_next - min_coord) > max_normal_isize_next:
                return True
            if max_coord_next < max_coord and (max_coord + alen - max_coord_next) > max_normal_isize_next:
                return True
            if max_coord_next >= max_coord and (max_coord_next + alen_next - max_coord) > max_normal_isize:
                return True
            return False
    return False
