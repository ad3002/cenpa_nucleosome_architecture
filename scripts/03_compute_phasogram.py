#!/usr/bin/env python3
"""
03_compute_phasogram.py

Computes spatial autocorrelation (phasogram) of nucleosome dyads along unbroken
centromeric alpha-satellite arrays to determine nucleosome repeat length (NRL)
inside and outside the Centromere Dip Region (CDR).

Uses an exact distance-bounded sliding window (no arbitrary neighbor truncation)
and supports parameterized fragment-length gates (default: 130-175 bp mononucleosome gate).
"""

import sys
import os
import math
import argparse
import subprocess
import collections

def compute_phasogram(bam_file, cdr_bed, out_tsv, min_len=130, max_len=175, max_lag=1200, bin_size=5):
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

    print(f"Streaming BAM pairs with gate [{min_len}, {max_len}] bp from {bam_file}...")
    cmd = ["samtools", "view", "-f", "2", "-F", "2304", bam_file]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1048576)

    dyads_cdr = collections.defaultdict(list)
    dyads_noncdr = collections.defaultdict(list)
    n_cdr_frags = 0
    n_noncdr_frags = 0

    for line in proc.stdout:
        p = line.rstrip("\n").split("\t")
        if len(p) < 9:
            continue
        tlen = int(p[8])
        if not (min_len <= tlen <= max_len):
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
                if cs <= genomic_mid < ce:
                    is_cdr = True
            if is_cdr:
                dyads_cdr[arr].append(dyad)
                n_cdr_frags += 1
            else:
                dyads_noncdr[arr].append(dyad)
                n_noncdr_frags += 1

    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError(f"samtools command failed with exit code {proc.returncode}")

    print(f"Collected {n_cdr_frags:,} CDR dyads and {n_noncdr_frags:,} Non-CDR dyads in gate [{min_len}, {max_len}] bp.")

    # Spatial Autocorrelation using exact distance-bounded window
    def autocorrelate_dyads(dyad_dict, max_dist, bsize):
        phas = collections.Counter()
        for arr, dlist in dyad_dict.items():
            if len(dlist) < 2:
                continue
            dlist.sort()
            n = len(dlist)
            for i in range(n):
                j = i + 1
                while j < n:
                    dist = dlist[j] - dlist[i]
                    if dist > max_dist:
                        break
                    bin_idx = int(math.floor((dist + bsize / 2.0) / float(bsize))) * bsize
                    phas[bin_idx] += 1
                    j += 1
        return phas

    print(f"Computing pairwise spatial autocorrelation up to {max_lag} bp...")
    phas_cdr = autocorrelate_dyads(dyads_cdr, max_lag, bin_size)
    phas_noncdr = autocorrelate_dyads(dyads_noncdr, max_lag, bin_size)

    # Write output TSV
    os.makedirs(os.path.dirname(os.path.abspath(out_tsv)), exist_ok=True)
    with open(out_tsv, "w") as out:
        out.write("distance_bp\tcdr_count\tnoncdr_count\n")
        for d in range(0, max_lag + bin_size, bin_size):
            out.write(f"{d}\t{phas_cdr[d]}\t{phas_noncdr[d]}\n")

    print(f"Phasogram written to {out_tsv} (range: 0 to {max_lag} bp, step {bin_size} bp)")

def main():
    parser = argparse.ArgumentParser(description="Compute spatial phasogram from paired-end BAM.")
    parser.add_argument("bam_file", help="Input coordinate-sorted BAM file")
    parser.add_argument("cdr_bed", help="BED file of CDR intervals")
    parser.add_argument("out_tsv", help="Output TSV path")
    parser.add_argument("--min-len", type=int, default=130, help="Minimum fragment length (default: 130)")
    parser.add_argument("--max-len", type=int, default=175, help="Maximum fragment length (default: 175)")
    parser.add_argument("--max-lag", type=int, default=1200, help="Maximum autocorrelation lag in bp (default: 1200)")
    parser.add_argument("--bin-size", type=int, default=5, help="Lag bin size in bp (default: 5)")
    args = parser.parse_args()

    compute_phasogram(
        args.bam_file,
        args.cdr_bed,
        args.out_tsv,
        min_len=args.min_len,
        max_len=args.max_len,
        max_lag=args.max_lag,
        bin_size=args.bin_size
    )

if __name__ == "__main__":
    main()
