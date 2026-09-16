#!/usr/bin/env python3
"""
03_compute_phasogram.py

Computes spatial autocorrelation (phasogram) of nucleosome dyads along unbroken
centromeric alpha-satellite arrays to determine nucleosome repeat length (NRL)
inside and outside the Centromere Dip Region (CDR).
"""

import sys
import os
import subprocess
import collections

def compute_phasogram(bam_file, cdr_bed, out_tsv):
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

    cmd = ["samtools", "view", "-f", "2", "-F", "2304", bam_file]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1048576)

    dyads_cdr = collections.defaultdict(list)
    dyads_noncdr = collections.defaultdict(list)

    for line in proc.stdout:
        p = line.rstrip("\n").split("\t")
        if len(p) < 9:
            continue
        tlen = int(p[8])
        if not (110 <= tlen <= 180):
            continue

        pos = int(p[3]) - 1
        arr = p[2]
        dyad = pos + tlen / 2.0

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
                dyads_cdr[arr].append(dyad)
            else:
                dyads_noncdr[arr].append(dyad)

    proc.wait()

    # Autocorrelation
    phas_cdr = collections.Counter()
    for arr, dlist in dyads_cdr.items():
        dlist.sort()
        for i in range(len(dlist)):
            for j in range(i + 1, min(i + 120, len(dlist))):
                dist = dlist[j] - dlist[i]
                if dist > 800:
                    break
                phas_cdr[int(round(dist / 5.0) * 5)] += 1

    phas_noncdr = collections.Counter()
    for arr, dlist in dyads_noncdr.items():
        dlist.sort()
        for i in range(len(dlist)):
            for j in range(i + 1, min(i + 120, len(dlist))):
                dist = dlist[j] - dlist[i]
                if dist > 800:
                    break
                phas_noncdr[int(round(dist / 5.0) * 5)] += 1

    with open(out_tsv, "w") as out:
        out.write("distance_bp\tcdr_count\tnoncdr_count\n")
        for d in range(0, 805, 5):
            out.write(f"{d}\t{phas_cdr[d]}\t{phas_noncdr[d]}\n")

    print(f"Phasogram written to {out_tsv}")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: 03_compute_phasogram.py <bam_file> <cdr_bed> <out_tsv>")
        sys.exit(1)
    compute_phasogram(sys.argv[1], sys.argv[2], sys.argv[3])
