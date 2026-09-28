#!/usr/bin/env python3
"""
experiments/exp03_dyad_vs_ends_geometry/analyze_dyad_vs_ends.py

Computational Experiment 3:
Resolving the Geometric Anchor of CENP-B Coupling:
Does the Bipartite Peak Pattern (+55 bp and +90-100 bp) Anchor to the Particle Dyad Center
or to Fragment Cleavage Termini?

Competing Hypotheses:
- H1: Fixed particle center (dyad). Slope beta(d ~ L) approx 0. Length variation arises
      from symmetric unpeeling/trimming of terminal arms (s ~ -0.5L, e ~ +0.5L).
- H2: Fixed cleavage terminus. One end is locked (e.g. against CENP-B barrier),
      yielding beta = -0.5 (fixed 5' end) or beta = +0.5 (fixed 3' end).
- H3: Monomer box switching. The two peaks reflect nearest-box assignment transitions
      between adjacent monomers rather than alternative conformations relative to the same box.

Evaluates both CHM13 Rep 2 (SRR13278683) and independent biological replicate Rep 1 (SRR13278684).
"""

import os
import sys
import csv
import bisect
import collections
import subprocess
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

EXP_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_OUT_DIR = os.path.join(EXP_DIR, "data")
FIG_OUT_DIR = os.path.join(EXP_DIR, "figures")
os.makedirs(DATA_OUT_DIR, exist_ok=True)
os.makedirs(FIG_OUT_DIR, exist_ok=True)

BASE_DIR = os.path.abspath(os.path.join(EXP_DIR, "..", ".."))
RAW_DIR = os.path.join(BASE_DIR, "raw_cache")
BOX_TSV = os.path.join(RAW_DIR, "chm13_cenpb_boxes_coords.tsv")

REPLICATES = {
    "CHM13_Rep2": os.path.join(RAW_DIR, "SRR13278683_slice.sorted.bam"),
    "CHM13_Rep1": os.path.join(RAW_DIR, "replication", "SRR13278684_slice.sorted.bam")
}

def load_cenpb_boxes():
    print(f"Loading CENP-B box coordinates from {BOX_TSV}...")
    boxes = collections.defaultdict(list)
    with open(BOX_TSV) as f:
        r = csv.DictReader(f, delimiter='\t')
        for row in r:
            arr = row['array_id']
            s = int(row['box_start'])
            e = int(row['box_end'])
            strand = row['strand']
            mid = (s + e) / 2.0
            boxes[arr].append((mid, s, e, strand))
    for arr in boxes:
        boxes[arr].sort(key=lambda x: x[0])
    total_boxes = sum(len(v) for v in boxes.values())
    print(f"Loaded {total_boxes:,} boxes across {len(boxes)} arrays.")
    return boxes

def extract_fragments(bam_path, boxes):
    print(f"Streaming primary proper pairs from {bam_path}...")
    cmd = ["samtools", "view", "-f", "2", "-F", "2304", bam_path]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1048576)
    
    records = []
    # record format:
    # (L, raw_dist, signed_dyad_offset, start_offset, end_offset, dist_box2, arr, pos)
    
    for line in proc.stdout:
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 9:
            continue
        tlen = int(parts[8])
        if tlen <= 0 or tlen < 70 or tlen > 200:
            continue
        pos = int(parts[3]) - 1
        arr = parts[2]
        if arr not in boxes or not boxes[arr]:
            continue
            
        L = tlen
        s = pos
        e = s + L
        m = s + L / 2.0
        
        arr_b = boxes[arr]
        idx = bisect.bisect_left(arr_b, (m,))
        cand_indices = []
        for offset in [-2, -1, 0, 1, 2]:
            j = idx + offset
            if 0 <= j < len(arr_b):
                cand_indices.append(j)
                
        if not cand_indices:
            continue
            
        sorted_cands = sorted(cand_indices, key=lambda j: abs(m - arr_b[j][0]))
        best_box = arr_b[sorted_cands[0]]
        b_mid, b_s, b_e, b_strand = best_box
        
        dist_box2 = abs(m - arr_b[sorted_cands[1]][0]) if len(sorted_cands) > 1 else np.nan
        
        raw_dist = abs(m - b_mid)
        # Orient relative to CENP-B box strand:
        if b_strand == '+':
            dyad_offset = m - b_mid      # positive: dyad downstream of box
            start_offset = s - b_mid     # fragment 5' end relative to box
            end_offset = e - b_mid       # fragment 3' end relative to box
        else:
            dyad_offset = b_mid - m      # positive: dyad downstream along reverse strand
            start_offset = b_mid - e     # 5' end along reverse strand is genomic e
            end_offset = b_mid - s       # 3' end along reverse strand is genomic s
            
        records.append({
            "length": L,
            "raw_dist": raw_dist,
            "dyad_offset": dyad_offset,
            "start_offset": start_offset,
            "end_offset": end_offset,
            "dist_box2": dist_box2,
            "array": arr,
            "pos": pos,
            "strand": b_strand
        })
        
    proc.stdout.close()
    proc.wait()
    print(f"Extracted {len(records):,} quality fragments.")
    return records

def bootstrap_regression(x, y, n_boot=500, ci=95):
    """Computes OLS slope and Theil-Sen slope with bootstrap confidence intervals."""
    n = len(x)
    ols_slope, ols_intercept = np.polyfit(x, y, 1)
    
    # Fast Theil-Sen robust estimator on representative subsample if n > 2000
    np.random.seed(42)
    if n > 2000:
        sub_idx = np.random.choice(n, size=2000, replace=False)
        res_ts = stats.theilslopes(y[sub_idx], x[sub_idx], alpha=0.95)
    else:
        res_ts = stats.theilslopes(y, x, alpha=0.95)
    ts_slope = res_ts[0]
    
    boot_slopes = []
    for _ in range(n_boot):
        b_idx = np.random.choice(n, size=n, replace=True)
        xb = x[b_idx]
        yb = y[b_idx]
        b_s, _ = np.polyfit(xb, yb, 1)
        boot_slopes.append(b_s)
        
    lower = np.percentile(boot_slopes, (100 - ci) / 2)
    upper = np.percentile(boot_slopes, 100 - (100 - ci) / 2)
    std_err = np.std(boot_slopes)
    
    return {
        "ols_slope": ols_slope,
        "ts_slope": ts_slope,
        "ci_lower": lower,
        "ci_upper": upper,
        "std_err": std_err
    }

def analyze_replicate(rep_name, records):
    print(f"\n--- Analyzing {rep_name} ---")
    
    # Convert to numpy arrays for vectorized masking
    lengths = np.array([r["length"] for r in records])
    raw_dists = np.array([r["raw_dist"] for r in records])
    dyad_offsets = np.array([r["dyad_offset"] for r in records])
    start_offsets = np.array([r["start_offset"] for r in records])
    end_offsets = np.array([r["end_offset"] for r in records])
    dist_box2 = np.array([r["dist_box2"] for r in records])
    
    # Windows:
    # Peak 1 (centered at ~55 bp, window [45, 65] bp)
    # Peak 2 (centered at ~95 bp, window [80, 110] bp)
    windows = {
        "Peak1_unsigned": (raw_dists >= 45) & (raw_dists <= 65),
        "Peak1_pos": (dyad_offsets >= 45) & (dyad_offsets <= 65),
        "Peak1_neg": (dyad_offsets >= -65) & (dyad_offsets <= -45),
        "Peak2_unsigned": (raw_dists >= 80) & (raw_dists <= 110),
        "Peak2_pos": (dyad_offsets >= 80) & (dyad_offsets <= 110),
        "Peak2_neg": (dyad_offsets >= -110) & (dyad_offsets <= -80),
    }
    
    results = {}
    for wname, mask in windows.items():
        n_frags = int(np.sum(mask))
        if n_frags < 100:
            print(f"Skipping {wname}: insufficient fragments ({n_frags})")
            continue
            
        sub_L = lengths[mask]
        sub_d = dyad_offsets[mask] if "pos" in wname or "neg" in wname else raw_dists[mask]
        sub_s = start_offsets[mask]
        sub_e = end_offsets[mask]
        
        dyad_stats = bootstrap_regression(sub_L, sub_d)
        start_stats = bootstrap_regression(sub_L, sub_s)
        end_stats = bootstrap_regression(sub_L, sub_e)
        
        # Test against H1 (slope = 0), H2_left (slope = -0.5), H2_right (slope = +0.5)
        # Compute z-scores:
        z_h1 = (dyad_stats["ols_slope"] - 0.0) / dyad_stats["std_err"]
        z_h2_left = (dyad_stats["ols_slope"] - (-0.5)) / dyad_stats["std_err"]
        z_h2_right = (dyad_stats["ols_slope"] - (+0.5)) / dyad_stats["std_err"]
        
        results[wname] = {
            "n_fragments": n_frags,
            "mean_length": float(np.mean(sub_L)),
            "std_length": float(np.std(sub_L)),
            "mean_dist": float(np.mean(sub_d)),
            "std_dist": float(np.std(sub_d)),
            "dyad_slope": dyad_stats["ols_slope"],
            "dyad_slope_ci95": [dyad_stats["ci_lower"], dyad_stats["ci_upper"]],
            "dyad_ts_slope": dyad_stats["ts_slope"],
            "start_slope": start_stats["ols_slope"],
            "end_slope": end_stats["ols_slope"],
            "z_score_vs_H1_center_fixed": z_h1,
            "z_score_vs_H2_5prime_fixed": z_h2_left,
            "z_score_vs_H2_3prime_fixed": z_h2_right,
            "p_val_vs_H1": 2 * (1 - stats.norm.cdf(abs(z_h1))),
            "p_val_vs_H2_5prime": 2 * (1 - stats.norm.cdf(abs(z_h2_left))),
            "p_val_vs_H2_3prime": 2 * (1 - stats.norm.cdf(abs(z_h2_right)))
        }
        
        print(f"[{wname}] N={n_frags:,}:")
        print(f"  Dyad slope beta = {dyad_stats['ols_slope']:+.4f} (95% CI: [{dyad_stats['ci_lower']:+.4f}, {dyad_stats['ci_upper']:+.4f}])")
        print(f"  Start slope beta = {start_stats['ols_slope']:+.4f}, End slope beta = {end_stats['ols_slope']:+.4f}")
        print(f"  Test vs H1 (beta=0): z={z_h1:.2f} | Test vs H2 (beta=-0.5): z={z_h2_left:.2f} | Test vs H2 (beta=+0.5): z={z_h2_right:.2f}")

    # Median by length curves for plotting
    median_curves = {}
    for wname in ["Peak1_unsigned", "Peak2_unsigned", "Peak1_pos", "Peak2_pos"]:
        if wname not in windows: continue
        mask = windows[wname]
        sub_L = lengths[mask]
        sub_d = dyad_offsets[mask] if "pos" in wname else raw_dists[mask]
        
        l_bins = range(100, 161, 5)
        curve_l = []
        curve_med = []
        curve_q25 = []
        curve_q75 = []
        for lb in l_bins:
            l_mask = (sub_L >= lb - 2) & (sub_L <= lb + 2)
            if np.sum(l_mask) >= 15:
                curve_l.append(lb)
                curve_med.append(float(np.median(sub_d[l_mask])))
                curve_q25.append(float(np.percentile(sub_d[l_mask], 25)))
                curve_q75.append(float(np.percentile(sub_d[l_mask], 75)))
        median_curves[wname] = {
            "lengths": curve_l, "medians": curve_med, "q25": curve_q25, "q75": curve_q75
        }
        
    return results, median_curves, records

def plot_experiment_figures(rep2_results, rep2_curves, rep2_records,
                            rep1_results, rep1_curves, rep1_records):
    print("\nGenerating Experiment 3 Publication Figures...")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    
    # Panel A: 2D Joint Density P(L, signed dyad offset) for Rep 2
    ax = axes[0, 0]
    r2_L = np.array([r["length"] for r in rep2_records])
    r2_d = np.array([r["dyad_offset"] for r in rep2_records])
    
    mask = (r2_L >= 90) & (r2_L <= 170) & (r2_d >= -150) & (r2_d <= 150)
    h = ax.hist2d(r2_L[mask], r2_d[mask], bins=[40, 60], cmap="magma_r", cmin=1)
    cbar = plt.colorbar(h[3], ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Fragment Count", fontsize=10)
    
    # Overlay Peak lines
    ax.axhline(55, color="#c2410c", linestyle="--", linewidth=1.5, label="Peak 1 (+55 bp gyre edge)")
    ax.axhline(-55, color="#c2410c", linestyle="--", linewidth=1.5)
    ax.axhline(95, color="#0369a1", linestyle="--", linewidth=1.5, label="Peak 2 (+95 bp linker)")
    ax.axhline(-95, color="#0369a1", linestyle="--", linewidth=1.5)
    ax.axhline(0, color="gray", linestyle=":", linewidth=1, label="Central Dyad (Depleted)")
    
    ax.set_title("A. Joint Distribution P(Length, Signed Dyad Offset) [CHM13 Rep 2]", fontsize=12, fontweight="bold")
    ax.set_xlabel("Fragment Length L (bp)", fontsize=11)
    ax.set_ylabel("Signed Dyad Offset from Box Center (bp)", fontsize=11)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax.grid(True, linestyle=":", alpha=0.5)

    # Panel B: Median Offset vs Length and Hypothesis Contrast
    ax = axes[0, 1]
    # Peak 1 (Rep 2)
    p1 = rep2_curves["Peak1_unsigned"]
    ax.errorbar(p1["lengths"], p1["medians"], 
                yerr=[np.array(p1["medians"]) - np.array(p1["q25"]), np.array(p1["q75"]) - np.array(p1["medians"])],
                fmt='o-', color="#c2410c", label="Peak 1 Observed (|d| ~ 55 bp)", linewidth=2, capsize=3)
    # Peak 2 (Rep 2)
    p2 = rep2_curves["Peak2_unsigned"]
    ax.errorbar(p2["lengths"], p2["medians"], 
                yerr=[np.array(p2["medians"]) - np.array(p2["q25"]), np.array(p2["q75"]) - np.array(p2["medians"])],
                fmt='s-', color="#0369a1", label="Peak 2 Observed (|d| ~ 95 bp)", linewidth=2, capsize=3)
    
    # Theoretical slopes
    ref_L = np.linspace(100, 160, 50)
    # H1: beta = 0
    ax.plot(ref_L, [55]*len(ref_L), 'k--', linewidth=1.5, label=r"H1 Model: Center Fixed ($\beta = 0$)")
    ax.plot(ref_L, [95]*len(ref_L), 'k--', linewidth=1.5)
    # H2: beta = -0.5 and +0.5
    ax.plot(ref_L, 55 - 0.5 * (ref_L - 130), 'r:', linewidth=1.5, label=r"H2 Model: 5' Fixed ($\beta = -0.5$)")
    ax.plot(ref_L, 95 + 0.5 * (ref_L - 130), 'b:', linewidth=1.5, label=r"H2 Model: 3' Fixed ($\beta = +0.5$)")
    
    ax.set_title("B. Geometric Anchor Test: Observed Slope vs H1 / H2", fontsize=12, fontweight="bold")
    ax.set_xlabel("Fragment Length L (bp)", fontsize=11)
    ax.set_ylabel("Distance to CENP-B Box Center (bp)", fontsize=11)
    ax.set_ylim(35, 125)
    ax.legend(loc="upper left", fontsize=8, framealpha=0.9)
    ax.grid(True, linestyle=":", alpha=0.5)

    # Panel C: Replicate Concordance (Rep 1 vs Rep 2 Slope Comparison)
    ax = axes[1, 0]
    categories = ["Peak 1\nUnsigned", "Peak 1\n(+55 bp)", "Peak 2\nUnsigned", "Peak 2\n(+95 bp)"]
    keys = ["Peak1_unsigned", "Peak1_pos", "Peak2_unsigned", "Peak2_pos"]
    
    rep2_slopes = [rep2_results[k]["dyad_slope"] for k in keys]
    rep2_errs = [rep2_results[k]["dyad_slope"] - rep2_results[k]["dyad_slope_ci95"][0] for k in keys]
    
    rep1_slopes = [rep1_results[k]["dyad_slope"] for k in keys]
    rep1_errs = [rep1_results[k]["dyad_slope"] - rep1_results[k]["dyad_slope_ci95"][0] for k in keys]
    
    x = np.arange(len(categories))
    width = 0.35
    
    rects1 = ax.bar(x - width/2, rep2_slopes, width, yerr=rep2_errs, label='CHM13 Rep 2 (Discovery)', 
                    color="#c2410c", alpha=0.85, capsize=4)
    rects2 = ax.bar(x + width/2, rep1_slopes, width, yerr=rep1_errs, label='CHM13 Rep 1 (Biological Replicate)', 
                    color="#f97316", alpha=0.85, capsize=4)
    
    ax.axhline(0, color="black", linestyle="-", linewidth=1.2)
    ax.axhline(-0.5, color="red", linestyle=":", linewidth=1.2, label=r"H2 Expectation ($\beta = -0.5$)")
    ax.axhline(+0.5, color="blue", linestyle=":", linewidth=1.2, label=r"H2 Expectation ($\beta = +0.5$)")
    
    ax.set_title("C. Cross-Replicate Validation of Slope Invariance", fontsize=12, fontweight="bold")
    ax.set_ylabel("Regression Slope beta (bp / bp of L)", fontsize=11)
    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=10)
    ax.set_ylim(-0.7, 0.7)
    ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
    ax.grid(True, linestyle=":", alpha=0.5)

    # Panel D: Test of Hypothesis H3 (Nearest Box vs 2nd Nearest Box Distance)
    ax = axes[1, 1]
    raw_d1 = np.array([r["raw_dist"] for r in rep2_records])
    raw_d2 = np.array([r["dist_box2"] for r in rep2_records if not np.isnan(r["dist_box2"])])
    
    # Filter to particles in Peak 1 vs Peak 2
    d2_p1 = [r["dist_box2"] for r in rep2_records if 45 <= r["raw_dist"] <= 65 and not np.isnan(r["dist_box2"])]
    d2_p2 = [r["dist_box2"] for r in rep2_records if 80 <= r["raw_dist"] <= 110 and not np.isnan(r["dist_box2"])]
    
    ax.hist(d2_p1, bins=range(0, 500, 10), density=True, alpha=0.6, color="#c2410c", label="Peak 1 particles: dist to 2nd box")
    ax.hist(d2_p2, bins=range(0, 500, 10), density=True, alpha=0.6, color="#0369a1", label="Peak 2 particles: dist to 2nd box")
    
    ax.axvline(340 - 55, color="#c2410c", linestyle="--", linewidth=1.5, label="Predicted 2nd box for Peak 1 (285 bp)")
    ax.axvline(340 - 95, color="#0369a1", linestyle="--", linewidth=1.5, label="Predicted 2nd box for Peak 2 (245 bp)")
    
    ax.set_title("D. Test of H3: Distance to 2nd Nearest CENP-B Box", fontsize=12, fontweight="bold")
    ax.set_xlabel("Distance to 2nd Nearest CENP-B Box (bp)", fontsize=11)
    ax.set_ylabel("Probability Density", fontsize=11)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax.grid(True, linestyle=":", alpha=0.5)
    
    plt.tight_layout()
    
    png_path = os.path.join(FIG_OUT_DIR, "Fig_Exp03_dyad_vs_ends_geometry.png")
    svg_path = os.path.join(FIG_OUT_DIR, "Fig_Exp03_dyad_vs_ends_geometry.svg")
    pdf_path = os.path.join(FIG_OUT_DIR, "Fig_Exp03_dyad_vs_ends_geometry.pdf")
    
    plt.savefig(png_path, dpi=300)
    plt.savefig(svg_path)
    plt.savefig(pdf_path)
    plt.close()
    print(f"Saved publication figures to {png_path}, {svg_path}, {pdf_path}")

def save_summary_tsvs(rep2_results, rep1_results):
    tsv_path = os.path.join(DATA_OUT_DIR, "exp03_slopes_and_hypothesis_tests.tsv")
    print(f"Saving statistical summary table to {tsv_path}...")
    
    fieldnames = [
        "cohort", "window", "n_fragments", "mean_length", "std_length",
        "mean_dist", "std_dist", "dyad_slope", "ci95_lower", "ci95_upper",
        "dyad_ts_slope", "start_slope", "end_slope", 
        "z_score_vs_H1_center", "p_val_vs_H1",
        "z_score_vs_H2_5prime", "p_val_vs_H2_5prime",
        "z_score_vs_H2_3prime", "p_val_vs_H2_3prime", "supported_hypothesis"
    ]
    
    with open(tsv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        
        for cname, res in [("CHM13_Rep2", rep2_results), ("CHM13_Rep1", rep1_results)]:
            for wname, row in res.items():
                # Hypothesis decision:
                # If slope is within [-0.25, 0.25] and rejects -0.5 and +0.5 with p < 1e-4 -> H1 (Center Fixed)
                # If slope is close to -0.5 -> H2 (5' Fixed)
                # If slope is close to +0.5 -> H2 (3' Fixed)
                s = row["dyad_slope"]
                if abs(s - 0.0) < 0.25 and row["p_val_vs_H2_5prime"] < 1e-4 and row["p_val_vs_H2_3prime"] < 1e-4:
                    decision = "H1 (Fixed Center / Symmetric Arms)"
                elif abs(s - (-0.5)) < 0.2:
                    decision = "H2 (5' Terminus Fixed)"
                elif abs(s - (+0.5)) < 0.2:
                    decision = "H2 (3' Terminus Fixed)"
                else:
                    decision = "Intermediate / Mixed"
                    
                w.writerow({
                    "cohort": cname,
                    "window": wname,
                    "n_fragments": row["n_fragments"],
                    "mean_length": f"{row['mean_length']:.2f}",
                    "std_length": f"{row['std_length']:.2f}",
                    "mean_dist": f"{row['mean_dist']:.2f}",
                    "std_dist": f"{row['std_dist']:.2f}",
                    "dyad_slope": f"{row['dyad_slope']:.4f}",
                    "ci95_lower": f"{row['dyad_slope_ci95'][0]:.4f}",
                    "ci95_upper": f"{row['dyad_slope_ci95'][1]:.4f}",
                    "dyad_ts_slope": f"{row['dyad_ts_slope']:.4f}",
                    "start_slope": f"{row['start_slope']:.4f}",
                    "end_slope": f"{row['end_slope']:.4f}",
                    "z_score_vs_H1_center": f"{row['z_score_vs_H1_center_fixed']:.2f}",
                    "p_val_vs_H1": f"{row['p_val_vs_H1']:.4e}",
                    "z_score_vs_H2_5prime": f"{row['z_score_vs_H2_5prime_fixed']:.2f}",
                    "p_val_vs_H2_5prime": f"{row['p_val_vs_H2_5prime']:.4e}",
                    "z_score_vs_H2_3prime": f"{row['z_score_vs_H2_3prime_fixed']:.2f}",
                    "p_val_vs_H2_3prime": f"{row['p_val_vs_H2_3prime']:.4e}",
                    "supported_hypothesis": decision
                })
                
    json_path = os.path.join(DATA_OUT_DIR, "exp03_results_summary.json")
    with open(json_path, "w") as f:
        json.dump({"CHM13_Rep2": rep2_results, "CHM13_Rep1": rep1_results}, f, indent=2)
    print(f"Saved machine-readable JSON to {json_path}.")

def main():
    boxes = load_cenpb_boxes()
    
    # 1. Process CHM13 Rep 2
    r2_records = extract_fragments(REPLICATES["CHM13_Rep2"], boxes)
    rep2_results, rep2_curves, r2_recs = analyze_replicate("CHM13_Rep2", r2_records)
    
    # 2. Process CHM13 Rep 1
    r1_records = extract_fragments(REPLICATES["CHM13_Rep1"], boxes)
    rep1_results, rep1_curves, r1_recs = analyze_replicate("CHM13_Rep1", r1_records)
    
    # 3. Generate figures & tables
    plot_experiment_figures(rep2_results, rep2_curves, r2_recs,
                            rep1_results, rep1_curves, r1_recs)
    save_summary_tsvs(rep2_results, rep1_results)
    
    print("\n=======================================================")
    print("  EXPERIMENT 3 EXECUTION COMPLETE!")
    print("=======================================================")

if __name__ == "__main__":
    main()
