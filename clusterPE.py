# Paired-ends Search_cliques
from utils import *
from utilsPE import *
import math


def search_cliquesPE(config: Config, comp_edge: TCompEdgeList, bam_record: List[PEAlignRecord],
                     svs: List[StructuralVariantRecord], svt: int):
    for comp_key, edge_list in comp_edge.items():
        edge_list.sort(key=cmp_to_key(sort_edge_records))
        it_w_edge = iter(edge_list)

        clique = set()
        incompatible = set()
        sv_start = -1
        sv_end = -1
        wiggle = 0
        cluster_ref_id = bam_record[next(it_w_edge).source].tid
        cluster_mate_ref_id = bam_record[next(it_w_edge).source].mtid
        sv_start, sv_end, wiggle = init_clique(bam_record[next(it_w_edge).source], sv_start=sv_start, sv_end=sv_end,
                                               wiggle=wiggle, svt=svt)

        if cluster_ref_id == cluster_mate_ref_id and sv_start >= sv_end:
            continue

        clique.add(next(it_w_edge).source)
        clique_grow = True

        while clique_grow:
            clique_grow = False
            for edge in edge_list:
                v = None
                if edge.source not in clique and edge.target in clique:
                    v = edge.source
                elif edge.source in clique and edge.target not in clique:
                    v = edge.target

                if v is None or v in incompatible:
                    continue

                clique_grow, sv_start, sv_end, wiggle = update_clique(bam_record[v], sv_start, sv_end, wiggle, svt)

                if clique_grow:
                    clique.add(v)
                else:
                    incompatible.add(v)
        # len(clique) = 91
        if len(clique) >= config.min_clique_size and sv_size_check(sv_start, sv_end, svt):

            sv_record = StructuralVariantRecord()
            sv_record.chr = cluster_ref_id
            sv_record.chr2 = cluster_mate_ref_id
            sv_record.sv_start = sv_start + 1
            sv_record.sv_end = sv_end + 1
            sv_record.pe_support = len(clique)

            ci_wiggle = max(abs(wiggle), 50)
            sv_record.ciposlow = -ci_wiggle
            sv_record.ciposhigh = ci_wiggle
            sv_record.ciendlow = -ci_wiggle
            sv_record.ciendhigh = ci_wiggle

            sv_record.mapq = 0
            mapQV = []
            for itC in clique:
                mapQV.append(bam_record[itC].map_quality)
                sv_record.mapq += bam_record[itC].map_quality
            mapQV.sort()
            sv_record.pe_map_quality = mapQV[len(mapQV) // 2]
            sv_record.sr_support = 0
            sv_record.sr_align_quality = 0
            sv_record.precise = False
            sv_record.svt = svt
            sv_record.insLen = 0
            sv_record.homLen = 0
            svs.append(sv_record)


def clusterPE(c: Config, bam_record: List[PEAlignRecord], svs: List[StructuralVariantRecord], varisize: int, svt: int):
    # Initialize variables
    comp = [0] * len(bam_record)
    num_comp = 0
    # Edge lists for each component
    comp_edge = {}  # {int:[EdgeRecord1,EdgeRecord2...]}

    # Initialize iterators and counters
    last_connected_node = 0
    last_connected_node_start = 0

    for i in range(len(bam_record)):
        if i > last_connected_node:
            if comp_edge:
                search_cliquesPE(c, comp_edge, bam_record, svs, svt)
                last_connected_node_start = last_connected_node
                comp_edge.clear()
        min_coord = min_coord_fn(bam_record[i].pos, bam_record[i].mpos, svt)
        max_coord = max_coord_fn(bam_record[i].pos, bam_record[i].mpos, svt)

        for j in range(i + 1, len(bam_record)):
            minCoord_next = min_coord_fn(bam_record[j].pos, bam_record[j].mpos, svt)
            if not (abs(minCoord_next + bam_record[j].alen - min_coord) <= varisize): break
            maxCoord_next = max_coord_fn(bam_record[j].pos, bam_record[j].mpos, svt)

            # Check that mate chromosomes agree (only for translocations)
            if bam_record[i].mtid != bam_record[j].mtid: continue
            # Check combinability of pairs
            if pairs_disagree(min_coord, max_coord, bam_record[i].alen, bam_record[i].max_normal_isize, minCoord_next, \
                              maxCoord_next, bam_record[j].alen, bam_record[j].max_normal_isize, svt):
                continue

            # Update last connected node
            if j > last_connected_node:
                last_connected_node = j

            # Assign components
            comp_index = 0
            if comp[i] == 0:
                if comp[j] == 0:
                    # Both vertices have no component
                    num_comp += 1
                    comp_index = num_comp
                    comp[i] = comp_index
                    comp[j] = comp_index
                    comp_edge[comp_index] = []
                else:
                    comp_index = comp[j]
                    comp[i] = comp_index
            else:
                if comp[j] == 0:
                    comp_index = comp[i]
                    comp[j] = comp_index
                else:
                    # Both vertices have a component
                    if comp[j] == comp[i]:
                        comp_index = comp[j]
                    else:
                        # Merge components
                        comp_index = comp[i]
                        other_index = comp[j]
                        if other_index < comp_index:
                            comp_index, other_index = other_index, comp_index
                        # Re-label other index
                        for k in range(last_connected_node_start, last_connected_node + 1):
                            if comp[k] == other_index:
                                comp[k] = comp_index
                        # Merge edge lists
                        comp_edge[comp_index].extend(comp_edge[other_index])
                        del comp_edge[other_index]

            # Append new edge
            if len(comp_edge[comp_index]) < c.graph_pruning:
                weight = int(math.log(abs(abs(minCoord_next - min_coord) - (maxCoord_next - max_coord)) - abs(
                    bam_record[i].median - bam_record[j].median) + 1) / math.log(2))
                comp_edge[comp_index].append(EdgeRecord(i, j, weight))
    if comp_edge:
        search_cliquesPE(c, comp_edge, bam_record, svs, svt)
        comp_edge.clear()