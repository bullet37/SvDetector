from utils import *
from utilsPE import *
import pysam
def hash_sr(rec):
    return hashlib.md5(str(rec).encode()).hexdigest()

def getVariability(config, lib):
    overall_variability = 0
    for library in lib:
        if library.maxNormalISize > overall_variability:
            overall_variability = library.maxNormalISize
        if library.rs > overall_variability:
            overall_variability = library.rs
    return overall_variability

def scanPE(c: Config, validRegions, svs: List[StructuralVariantRecord],
                srSVs: List[StructuralVariantRecord], srStore, sampleLib: LibraryInfo):
    samfile = [pysam.AlignmentFile(file, "rb") for file in c.files]
    hdr = samfile[0].header

    # Bam alignment record vector
    bamRecords = [ list() for _ in range(2 * DELLY_SVT_TRANS)]
    # Iterate all samples
    for file_c in range(len(c.files)):
        matetra = [{} for _ in range(len(c.files))]
        n_targets = samfile[file_c].nreferences
        # Iterate all chromosome for that sample
        for refIndex in range(n_targets):
            if not validRegions[refIndex]:
                continue
            nodata = True
            suffix = "cram"
            str1 = c.files[file_c]

            if str1.endswith(suffix) and str1[-len(suffix):] == suffix:
                nodata = False

            stats = samfile[file_c].get_index_statistics()
            mapped = stats[refIndex].mapped
            unmapped = stats[refIndex].unmapped
            if mapped:
                nodata = False
            if nodata:
                continue

            mateMap = {}
            # Read alignments
            for vRIt in validRegions[refIndex]:
                lastAlignedPosReads = set()
                lastAlignedPos = 0
                refName = samfile[file_c].references[refIndex]
                iter = samfile[file_c].fetch(refName, vRIt.lower(), vRIt.upper())
                for rec in iter:

                    if rec.flag & (pysam.FQCFAIL | pysam.FDUP | pysam.FUNMAP): continue
                    if rec.mapping_quality < c.minMapQual or rec.tid < 0: continue
                    # Paired-end clustering, build up bamRecord
                    seed = hash_sr(rec)
                    if rec.flag & (pysam.FPAIRED):

                        if sampleLib[file_c].median == 0: continue  # Single-end library
                        if rec.flag & (pysam.FSECONDARY | pysam.FSUPPLEMENTARY): continue  # Secondary/supplementary alignments
                        if rec.next_reference_id < 0 or rec.flag & pysam.FMUNMAP: continue  # Mate unmapped or blacklisted chr
                        if not validRegions[rec.next_reference_id]: continue
                        if is_rec_translocation(rec) and rec.mapping_quality < c.minTraQual:  continue

                        # SV type
                        svt = isizeMappingPos(rec, sampleLib[file_c].maxISizeCutoff)
                        if svt == -1: continue
                        if c.svtset and svt not in c.svtset: continue

                        # Check library-specific insert size for deletions
                        if svt == 2 and sampleLib[file_c].maxISizeCutoff > np.abs(rec.isize): continue

                        if rec.pos > lastAlignedPos:
                            lastAlignedPosReads.clear()
                            lastAlignedPos = rec.pos

                        if firstPairObs(rec, lastAlignedPosReads):
                            # First read
                            lastAlignedPosReads.add(seed)
                            hv = hash_pair(rec)
                            if is_translocation(svt):
                                matetra[file_c][hv] = (rec.mapping_quality, alignmentLength(rec))
                            else:
                                mateMap[hv] = (rec.mapping_quality, alignmentLength(rec))
                        else:
                          # Second read
                          hv = hash_pair_mate(rec)
                          alenmate = 0
                          pairQuality = 0
                          if is_translocation(svt):
                              if hv in matetra[file_c] and matetra[file_c][hv][0]:
                                  pairQuality = min(matetra[file_c][hv][0], rec.mapping_quality)
                                  alenmate = matetra[file_c][hv][1]
                                  matetra[file_c][hv] = (0, alenmate)
                          else:
                              if hv in mateMap and mateMap[hv][0]:
                                  pairQuality = min(mateMap[hv][0], rec.mapping_quality)
                                  alenmate = mateMap[hv][1]
                                  mateMap[hv] = (0, alenmate)

                          bamRecord = BamAlignRecord(rec, pairQuality, alignmentLength(rec), alenmate, sampleLib[file_c].median,
                                         sampleLib[file_c].mad, sampleLib[file_c].maxNormalISize)
                          bamRecords[svt].append(bamRecord)
                          sampleLib[file_c].abnormal_pairs += 1

    # Maximum variability in insert size
    varisize = getVariability(c, sampleLib)
    for svt in range(len(bamRecords)):
        if c.svtset and svt not in c.svtset: continue
        if not bamRecords[svt]: continue
        # Sort BAM records according to position
        bamRecords[svt].sort(key=cmp_to_key(SortBamRecords))
        # Cluster
        clusterPR(c, bamRecords[svt], svs, varisize, svt)

    return bamRecords,svs

