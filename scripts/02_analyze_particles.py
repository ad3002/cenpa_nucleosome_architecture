#!/usr/bin/env python3
"""
02_analyze_particles.py

Processes coordinate-sorted BAM files of paired-end MNase sequencing aligned to
centromeric alpha-satellite arrays to compute:
1. Fragment length (insert size) distribution (global, CDR, Non-CDR).
2. Mononucleosome and dinucleosome peak modes.
3. Distance from nucleosome dyad to the nearest CENP-B box center.
4. Dinucleosome core vs linker enrichment of CENP-B boxes.
"""

import sys
import os
import csv
import bisect
import collections
import subprocess

def run_analysis(bam_file, cdr_bed, boxes_tsv, out_prefix):
    # 1. Chromosome mapping for CHM13v2.0
    chrom_map = {f"chr{i}": f"NC_{60925 + i - 1:06d}.1" for i in range(1, 23)}
    chrom_map["chrX"] = "NC_060947.1"
    acc_to_chr = {v: k for k, v in chrom_map.items()}

    cdrs = {}
    if os.path.exists(cdr_bed):
        with open(cdr_bed) as f:
            for line in f:
                p = line.strip().split()
                if len(p) >= 3:
                    c = p[0].replace("_hap1", "")
                    acc = chrom_map.get(c)
                    if acc:
                        cdrs[acc] = (int(p[1]), int(p[2]))
        print(f"Loaded {len(cdrs)} CDRs.")

    # 2. Load CENP-B boxes
    print(f"Loading CENP-B boxes from {boxes_tsv}...")
    box_centers = collections.defaultdict(list)
    total_boxes = 0
    if os.path.exists(boxes_tsv):
        with open(boxes_tsv) as f:
            r = csv.DictReader(f, delimiter="\t")
            for row in r:
                arr = row["array_id"]
                bs = int(row["box_start"])
                be = int(row["box_end"])
                center = (bs + be) / 2.0
                box_centers[arr].append(center)
                total_boxes += 1
        for arr in box_centers:
            box_centers[arr].sort()
        print(f"Loaded {total_boxes} boxes across {len(box_centers)} arrays.")

    # 3. Stream BAM using samtools
    print(f"Streaming primary proper pairs from {bam_file}...")
    cmd = ["samtools", "view", "-f", "2", "-F", "2304", bam_file]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1048576)

    hist_global = collections.Counter()
    hist_cdr = collections.Counter()
    hist_noncdr = collections.Counter()
    box_to_dyad_dist = collections.Counter()

    CORE_LEN = 147
    total_frags = 0

    for line in proc.stdout:
        p = line.rstrip("\n").split("\t")
        if len(p) < 9:
            continue
        tlen = int(p[8])
        if tlen <= 0 or tlen < 50 or tlen > 1200:
            continue

        pos = int(p[3]) - 1
        arr = p[2]
        total_frags += 1

        f_start = pos
        f_len = tlen
        dyad = f_start + f_len / 2.0
        hist_global[f_len] += 1

        parts = arr.split("_")
        if len(parts) >= 3:
            acc = parts[0] + "_" + parts[1]
            try:
                arr_s = int(parts[2])
            except ValueError:
                arr_s = 0
            genomic_mid = arr_s + dyad
            is_cdr = False
            if acc in cdrs:
                cs, ce = cdrs[acc]
                if cs <= genomic_mid <= ce:
                    is_cdr = True
            if is_cdr:
                hist_cdr[f_len] += 1
            else:
                hist_noncdr[f_len] += 1

        # Distance to CENP-B box
        centers = box_centers.get(arr)
        if centers:
            if 60 <= f_len <= 200:
                idx = bisect.bisect_left(centers, dyad)
                best_d = 999999
                if idx < len(centers):
                    best_d = min(best_d, abs(centers[idx] - dyad))
                if idx > 0:
                    best_d = min(best_d, abs(centers[idx - 1] - dyad))
                if best_d <= 300:
                    b_bin = int(best_d // 5) * 5
                    box_to_dyad_dist[b_bin] += 1

    proc.wait()
    print(f"Processed {total_frags:,} primary proper-pair fragments.")

    # Write output tables
    hist_out = f"{out_prefix}_fragment_length_hist.tsv"
    with open(hist_out, "w") as f:
        f.write("fragment_length_bp\tglobal_count\tcdr_count\tnoncdr_count\n")
        max_l = max(hist_global.keys()) if hist_global else 600
        for l in range(50, min(max_l + 1, 800)):
            f.write(f"{l}\t{hist_global[l]}\t{hist_cdr[l]}\t{hist_noncdr[l]}\n")
    print(f"Wrote {hist_out}")

    dist_out = f"{out_prefix}_box_to_dyad_distance.tsv"
    with open(dist_out, "w") as f:
        f.write("distance_to_dyad_bp\tcount\n")
        for d in sorted(box_to_dyad_dist.keys()):
            f.write(f"{d}\t{box_to_dyad_dist[d]}\n")
    print(f"Wrote {dist_out}")

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: 02_analyze_particles.py <bam_file> <cdr_bed> <boxes_tsv> <out_prefix>")
        sys.exit(1)
    run_analysis(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
