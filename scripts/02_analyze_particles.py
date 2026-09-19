#!/usr/bin/env python3
"""
02_analyze_particles.py

Processes coordinate-sorted BAM files of paired-end MNase sequencing aligned to
centromeric alpha-satellite arrays to compute:
1. Fragment length (insert size) distribution (global, CDR, Non-CDR).
2. Distance from nucleosome dyad to the nearest CENP-B box center.
3. Optional per-chromosome summary table (--summary-out).
"""

import sys
import os
import math
import csv
import bisect
import collections
import subprocess
import argparse

def run_analysis(bam_file, cdr_bed, boxes_tsv, hist_out, dist_out, summary_out=None, max_len=800):
    # 1. Chromosome mapping for CHM13v2.0
    chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
    chrom_map = {f"chr{i}": f"NC_{60925 + i - 1:06d}.1" for i in range(1, 23)}
    chrom_map["chrX"] = "NC_060947.1"
    acc_to_chr = {v: k for k, v in chrom_map.items()}

    def find_5bp_mode(counter, min_len=100, max_len=200):
        bin_counts = collections.Counter()
        for l, cnt in counter.items():
            if min_len <= l <= max_len:
                b = int(math.floor((l + 2.5) / 5.0)) * 5
                bin_counts[b] += cnt
        if not bin_counts or sum(bin_counts.values()) == 0:
            return "NA"
        return max(bin_counts, key=bin_counts.get)

    def find_di_peak_func(counter, min_len=250, max_len=350, min_support=50, snr_threshold=1.8, min_peak_count=15):
        bin_counts = collections.Counter()
        for l, cnt in counter.items():
            if min_len <= l <= max_len:
                b = int(math.floor((l + 2.5) / 5.0)) * 5
                bin_counts[b] += cnt
        if not bin_counts or sum(bin_counts.values()) == 0:
            return None
        tot_in_window = sum(bin_counts.values())
        if tot_in_window < min_support:
            return None
        possible_bins = list(range(min_len, max_len + 5, 5))
        avg_background = tot_in_window / float(len(possible_bins))
        top_b, top_c = bin_counts.most_common(1)[0]
        if top_c >= snr_threshold * avg_background and top_c >= min_peak_count:
            return top_b
        return None

    if not os.path.exists(bam_file):
        if os.path.exists(hist_out):
            print(f"Notice: BAM file {bam_file} not found; analyzing pre-computed histogram from {hist_out}...")
            hist_global = collections.Counter()
            hist_cdr = collections.Counter()
            hist_noncdr = collections.Counter()
            with open(hist_out) as f:
                for r in csv.DictReader(f, delimiter="\t"):
                    l = int(r["fragment_length_bp"])
                    hist_global[l] = int(r["global_count"])
                    hist_cdr[l] = int(r["cdr_count"])
                    hist_noncdr[l] = int(r["noncdr_count"])
            
            if summary_out and os.path.exists(summary_out):
                with open(summary_out) as f:
                    rows = list(csv.DictReader(f, delimiter="\t"))
                
                glob_m_c = find_5bp_mode(hist_cdr)
                glob_m_nc = find_5bp_mode(hist_noncdr)
                glob_di_c = find_di_peak_func(hist_cdr)
                glob_di_nc = find_di_peak_func(hist_noncdr)
                glob_nrl_c = (glob_di_c - glob_m_c) if (glob_di_c is not None and glob_m_c != "NA") else "None"
                glob_nrl_nc = (glob_di_nc - glob_m_nc) if (glob_di_nc is not None and glob_m_nc != "NA") else "None"
                glob_delta = (glob_nrl_c - glob_nrl_nc) if (isinstance(glob_nrl_c, int) and isinstance(glob_nrl_nc, int)) else "NA"
                
                for r in rows:
                    if r["chrom"] == "GLOBAL":
                        r["Mono_CDR"] = str(glob_m_c)
                        r["Di_CDR"] = str(glob_di_c if glob_di_c is not None else "None")
                        r["NRL_CDR"] = str(glob_nrl_c)
                        r["Mono_NonCDR"] = str(glob_m_nc)
                        r["Di_NonCDR"] = str(glob_di_nc if glob_di_nc is not None else "None")
                        r["NRL_NonCDR"] = str(glob_nrl_nc)
                        r["Delta_NRL"] = str(glob_delta)
                    else:
                        # Ensure Di and NRL match threshold if Di was called on tiny count
                        try:
                            # If Di was None, keep None
                            if r.get("Di_CDR") in ["None", "—", ""]:
                                r["Di_CDR"] = "None"
                                r["NRL_CDR"] = "None"
                            if r.get("Di_NonCDR") in ["None", "—", ""]:
                                r["Di_NonCDR"] = "None"
                                r["NRL_NonCDR"] = "None"
                            nc = int(r["NRL_CDR"])
                            nnc = int(r["NRL_NonCDR"])
                            r["Delta_NRL"] = str(nc - nnc)
                        except (ValueError, TypeError):
                            r["Delta_NRL"] = "NA"
                
                with open(summary_out, "w") as f:
                    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
                    writer.writeheader()
                    for r in rows:
                        writer.writerow(r)
                print(f"Updated {summary_out}")
            return
        raise FileNotFoundError(f"Input BAM file not found: {bam_file}")

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
        print(f"Loaded {len(cdrs)} CDR intervals.")
    else:
        print(f"Notice: CDR BED file {cdr_bed} not found; CDR partitioning will be skipped.")

    # 2. Load CENP-B boxes
    box_centers = collections.defaultdict(list)
    total_boxes = 0
    if os.path.exists(boxes_tsv):
        print(f"Loading CENP-B boxes from {boxes_tsv}...")
        with open(boxes_tsv) as f:
            r = csv.DictReader(f, delimiter="\t")
            for row in r:
                arr = row.get("array_id", row.get("chrom", ""))
                bs = int(row.get("box_start", row.get("start", 0)))
                be = int(row.get("box_end", row.get("end", 0)))
                center = (bs + be) / 2.0
                box_centers[arr].append(center)
                total_boxes += 1
        for arr in box_centers:
            box_centers[arr].sort()
        print(f"Loaded {total_boxes} boxes across {len(box_centers)} arrays.")
    else:
        print(f"Notice: CENP-B box coordinates file {boxes_tsv} not found; box distances will be empty.")

    # 3. Stream BAM using samtools
    print(f"Streaming primary proper pairs from {bam_file}...")
    cmd = ["samtools", "view", "-f", "2", "-F", "2304", bam_file]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1048576)

    hist_global = collections.Counter()
    hist_cdr = collections.Counter()
    hist_noncdr = collections.Counter()
    box_to_dyad_dist = collections.Counter()

    # Per-chromosome tracking
    chr_cdr_counts = collections.Counter()
    chr_noncdr_counts = collections.Counter()
    chr_cdr_hist = collections.defaultdict(collections.Counter)
    chr_noncdr_hist = collections.defaultdict(collections.Counter)

    total_frags = 0

    for line in proc.stdout:
        p = line.rstrip("\n").split("\t")
        if len(p) < 9:
            continue
        tlen = int(p[8])
        if tlen <= 0 or tlen < 50 or tlen > max_len:
            continue

        pos = int(p[3]) - 1
        arr = p[2]
        total_frags += 1

        f_start = pos
        f_len = tlen
        dyad = f_start + f_len / 2.0
        hist_global[f_len] += 1

        parts = arr.split("_")
        chr_name = None
        if len(parts) >= 3:
            acc = parts[0] + "_" + parts[1]
            chr_name = acc_to_chr.get(acc)
            try:
                arr_s = int(parts[2])
            except ValueError:
                arr_s = 0
            genomic_mid = arr_s + dyad
            is_cdr = False
            if acc in cdrs:
                cs, ce = cdrs[acc]
                # Half-open 0-based interval [start, end)
                if cs <= genomic_mid < ce:
                    is_cdr = True
            if is_cdr:
                hist_cdr[f_len] += 1
                if chr_name:
                    chr_cdr_counts[chr_name] += 1
                    chr_cdr_hist[chr_name][f_len] += 1
            else:
                hist_noncdr[f_len] += 1
                if chr_name:
                    chr_noncdr_counts[chr_name] += 1
                    chr_noncdr_hist[chr_name][f_len] += 1
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
    if proc.returncode != 0:
        raise RuntimeError(f"samtools execution failed with exit code {proc.returncode}")

    print(f"Processed {total_frags:,} primary proper-pair fragments.")

    # Write output tables
    os.makedirs(os.path.dirname(os.path.abspath(hist_out)), exist_ok=True)
    with open(hist_out, "w") as f:
        f.write("fragment_length_bp\tglobal_count\tcdr_count\tnoncdr_count\n")
        for l in range(50, max_len + 1):
            f.write(f"{l}\t{hist_global[l]}\t{hist_cdr[l]}\t{hist_noncdr[l]}\n")
    print(f"Wrote {hist_out}")

    os.makedirs(os.path.dirname(os.path.abspath(dist_out)), exist_ok=True)
    with open(dist_out, "w") as f:
        f.write("distance_to_dyad_bp\tcount\n")
        for d in range(0, 305, 5):
            f.write(f"{d}\t{box_to_dyad_dist[d]}\n")
    print(f"Wrote {dist_out}")

    if summary_out:
        def find_5bp_mode(counter, min_len=100, max_len=200):
            bin_counts = collections.Counter()
            for l, cnt in counter.items():
                if min_len <= l <= max_len:
                    b = int(round(l / 5.0)) * 5
                    bin_counts[b] += cnt
            if not bin_counts:
                return "NA"
            return max(bin_counts, key=bin_counts.get)

        def find_di_peak(counter, min_len=250, max_len=350, min_support=50, snr_threshold=1.8, min_peak_count=15):
            bin_counts = collections.Counter()
            for l, cnt in counter.items():
                if min_len <= l <= max_len:
                    b = int(round(l / 5.0)) * 5
                    bin_counts[b] += cnt
            if not bin_counts:
                return None
            tot_in_window = sum(bin_counts.values())
            # Require minimum statistical support in the dinucleosome search window
            if tot_in_window < min_support:
                return None
            # Evaluate average background across all possible 5-bp bins in the window
            possible_bins = list(range(min_len, max_len + 5, 5))
            avg_background = tot_in_window / float(len(possible_bins))
            top_b, top_c = bin_counts.most_common(1)[0]
            if top_c >= snr_threshold * avg_background and top_c >= min_peak_count:
                return top_b
            return None

        os.makedirs(os.path.dirname(os.path.abspath(summary_out)), exist_ok=True)
        with open(summary_out, "w") as f:
            f.write("chrom\tN_cdr\tMono_CDR\tDi_CDR\tNRL_CDR\tN_noncdr\tMono_NonCDR\tDi_NonCDR\tNRL_NonCDR\tDelta_NRL\n")
            for c in chrom_order:
                n_c = chr_cdr_counts[c]
                n_nc = chr_noncdr_counts[c]
                m_c = find_5bp_mode(chr_cdr_hist[c])
                m_nc = find_5bp_mode(chr_noncdr_hist[c])
                di_c = find_di_peak(chr_cdr_hist[c])
                di_nc = find_di_peak(chr_noncdr_hist[c])
                nrl_c = di_c - m_c if (di_c is not None and m_c != "NA") else "None"
                nrl_nc = di_nc - m_nc if (di_nc is not None and m_nc != "NA") else "None"
                delta = (nrl_c - nrl_nc) if (isinstance(nrl_c, int) and isinstance(nrl_nc, int)) else "NA"
                f.write(f"{c}\t{n_c}\t{m_c}\t{di_c if di_c is not None else 'None'}\t{nrl_c}\t{n_nc}\t{m_nc}\t{di_nc if di_nc is not None else 'None'}\t{nrl_nc}\t{delta}\n")
            
            # Global row
            tot_cdr = sum(chr_cdr_counts.values())
            tot_noncdr = sum(hist_noncdr.values())
            glob_m_c = find_5bp_mode(hist_cdr)
            glob_m_nc = find_5bp_mode(hist_noncdr)
            glob_di_c = find_di_peak(hist_cdr)
            glob_di_nc = find_di_peak(hist_noncdr)
            glob_nrl_c = (glob_di_c - glob_m_c) if (glob_di_c is not None and glob_m_c != "NA") else "None"
            glob_nrl_nc = (glob_di_nc - glob_m_nc) if (glob_di_nc is not None and glob_m_nc != "NA") else "None"
            glob_delta = (glob_nrl_c - glob_nrl_nc) if (isinstance(glob_nrl_c, int) and isinstance(glob_nrl_nc, int)) else "NA"
            f.write(f"GLOBAL\t{tot_cdr}\t{glob_m_c}\t{glob_di_c if glob_di_c is not None else 'None'}\t{glob_nrl_c}\t{tot_noncdr}\t{glob_m_nc}\t{glob_di_nc if glob_di_nc is not None else 'None'}\t{glob_nrl_nc}\t{glob_delta}\n")
        print(f"Wrote {summary_out}")

def main():
    parser = argparse.ArgumentParser(description="Analyze nucleosome particle sizes and CENP-B box distances from BAM.")
    parser.add_argument("bam_file", help="Input sorted BAM file")
    parser.add_argument("cdr_bed", help="BED file with CDR coordinates")
    parser.add_argument("boxes_tsv", help="TSV file with CENP-B box coordinates")
    parser.add_argument("out_prefix_or_hist", help="Output prefix or direct path to fragment length histogram TSV")
    parser.add_argument("dist_out", nargs="?", default=None, help="Optional direct path to box distance TSV")
    parser.add_argument("--summary-out", default=None, help="Optional output path for per-chromosome summary TSV")
    parser.add_argument("--max-len", type=int, default=800, help="Maximum fragment length to process (default: 800)")
    args = parser.parse_args()

    if args.dist_out:
        hist_out = args.out_prefix_or_hist
        dist_out = args.dist_out
    else:
        prefix = args.out_prefix_or_hist
        hist_out = f"{prefix}_fragment_length_hist.tsv"
        dist_out = f"{prefix}_box_to_dyad_distance.tsv"

    run_analysis(args.bam_file, args.cdr_bed, args.boxes_tsv, hist_out, dist_out, summary_out=args.summary_out, max_len=args.max_len)

if __name__ == "__main__":
    main()
