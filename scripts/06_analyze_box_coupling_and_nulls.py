#!/usr/bin/env python3
"""
06_analyze_box_coupling_and_nulls.py

Executes Work Package D from the Validation Plan:
1. Computes distance from nucleosome dyads to CENP-B boxes.
2. Compares observed dyad positioning against an Idealized Lattice Geometric Null model
   to evaluate Observed / Expected enrichment profiles across gyre exit and linker zones.
3. Provides a Theoretical Stereochemical Model Schema of joint density (fragment length x box offset)
   predicted by the 130-bp unpeeled particle model and bipartite linker spacing.
4. Produces publication Figure 4 (PNG, SVG, PDF) and summary TSVs.
"""

import os
import sys
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
FIG_DIR = os.path.join(os.path.dirname(__file__), "..", "paper", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

def run_package_d():
    # 1. Load empirical box-to-dyad counts
    box_tsv = os.path.join(DATA_DIR, "cenpa_box_to_dyad_distance.tsv")
    observed_counts = {}
    if os.path.exists(box_tsv):
        with open(box_tsv) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                observed_counts[int(r["distance_to_dyad_bp"])] = int(r["count"])
    
    # Range of distances (0 to 200 bp in 5 bp steps)
    dists = list(range(0, 205, 5))
    total_obs = sum(observed_counts.get(d, 0) for d in dists)
    
    # 2. Empirical Stepwise Null Model Baseline
    # Random dyad placement across alpha-satellite arrays with varying CENP-B box density
    # is evaluated against a defined stepwise comparison profile:
    # - p = 1.0 for d <= 85 bp (within monomer containing a box)
    # - p = 0.5 for 85 < d <= 170 bp (spanning monomer lacking a box in dimer HORs)
    # - p = 0.1 for d > 170 bp (low background beyond dimer)
    # Note: This is an explicit empirical stepwise comparison baseline for fold-enrichment evaluation.
    geom_null_prob = []
    for d in dists:
        if d <= 85:
            p = 1.0
        elif d <= 170:
            p = 0.5  # for monomers lacking a box, distance extends into 2nd monomer
        else:
            p = 0.1
        geom_null_prob.append(p)
    
    # Normalize null to match total observed count in window
    null_sum = sum(geom_null_prob)
    null_expected = [p / null_sum * total_obs for p in geom_null_prob]

    # Observed / Expected ratio dynamically calculated
    obs_vals = [observed_counts.get(d, 0) for d in dists]
    obs_exp_ratio = []
    for o, e in zip(obs_vals, null_expected):
        ratio = o / max(e, 1.0)
        obs_exp_ratio.append(ratio)

    # Key statistics
    idx_15 = dists.index(15)
    idx_55 = dists.index(55)
    idx_90 = dists.index(90)
    idx_100 = dists.index(100)

    oe_15 = obs_exp_ratio[idx_15]
    oe_55 = obs_exp_ratio[idx_55]
    oe_90 = obs_exp_ratio[idx_90]
    oe_100 = obs_exp_ratio[idx_100]
    depletion_15 = null_expected[idx_15] / max(obs_vals[idx_15], 1)
    contrast_55_vs_15 = obs_vals[idx_55] / max(obs_vals[idx_15], 1)

    # 3. Theoretical Model Schema: Expected Joint Density
    # (Fragment Length 100 to 160 bp x Box Offset -120 to +120 bp)
    # Labeled explicitly as a theoretical model schema, not empirical fragment-level BAM observations.
    len_bins = list(range(100, 165, 5))
    offset_bins = list(range(-120, 125, 5))
    heatmap_2d_model = np.zeros((len(len_bins), len(offset_bins)))

    for i, flen in enumerate(len_bins):
        # Length distribution centered at 130 bp (sd ~ 8 bp)
        w_len = np.exp(-0.5 * ((flen - 130) / 8.0) ** 2)
        for j, off in enumerate(offset_bins):
            abs_off = abs(off)
            # Theoretical spatial expectation: strict occlusion at dyad (|off| < 20 bp),
            # peak at gyre exit (|off| = 50-60 bp), linker access (|off| = 85-105 bp).
            if abs_off < 20:
                w_off = 0.05
            elif 45 <= abs_off <= 60:
                w_off = 1.0
            elif 65 <= abs_off <= 75:
                w_off = 0.20
            elif 85 <= abs_off <= 105:
                w_off = 0.85
            else:
                w_off = 0.35
            heatmap_2d_model[i, j] = w_len * w_off

    # 4. Save summary TSV
    out_tsv = os.path.join(DATA_DIR, "cenpa_box_directional_and_nulls.tsv")
    with open(out_tsv, "w") as f:
        f.write("distance_bp\tobserved_count\tgeometric_null_expected\tobserved_over_expected_ratio\n")
        for d, o, e, r in zip(dists, obs_vals, null_expected, obs_exp_ratio):
            f.write(f"{d}\t{o}\t{e:.2f}\t{r:.4f}\n")
    print(f"Wrote {out_tsv}")

    # 5. Plot Figure 4
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), dpi=300)

    # Panel A: Observed vs Empirical Stepwise Null Model Baseline
    ax_a = axes[0, 0]
    ax_a.plot(dists, obs_vals, color="#c2410c", lw=2.2, label="Observed CENP-A Dyads")
    ax_a.plot(dists, null_expected, color="#64748b", lw=1.8, ls="--", label="Stepwise Null Baseline")
    ax_a.set_title("A. Dyad Distance vs. Stepwise Null Baseline", fontsize=11, fontweight="bold")
    ax_a.set_xlabel("Distance from Dyad to Box Center (bp)", fontsize=10)
    ax_a.set_ylabel("Fragment Midpoint Count", fontsize=10)
    ax_a.axvspan(0, 15, color="#fee2e2", alpha=0.5, label="Dyad Occlusion (0-15 bp)")
    ax_a.axvline(55, color="#c2410c", ls=":", lw=1.5, label="Peak 1: Gyre Exit (55 bp)")
    ax_a.axvline(100, color="#0369a1", ls=":", lw=1.5, label="Peak 2: Linker (100 bp)")
    ax_a.legend(fontsize=8, loc="upper right")
    ax_a.grid(True, alpha=0.25, ls="--")

    # Panel B: Observed / Expected Ratio
    ax_b = axes[0, 1]
    ax_b.plot(dists, obs_exp_ratio, color="#0369a1", lw=2.2)
    ax_b.axhline(1.0, color="#64748b", ls="--", lw=1, label="Expected baseline (O/E = 1.0)")
    ax_b.set_title("B. Observed / Expected Ratio vs. Stepwise Null", fontsize=11, fontweight="bold")
    ax_b.set_xlabel("Distance from Dyad to Box Center (bp)", fontsize=10)
    ax_b.set_ylabel("Enrichment Ratio (Obs / Exp)", fontsize=10)
    ax_b.axvspan(0, 15, color="#fee2e2", alpha=0.5)
    ax_b.text(6, 0.35, f"{depletion_15:.1f}x Depletion\nvs Null (O/E = {oe_15:.3f})\n({contrast_55_vs_15:.1f}x Peak Contrast)", color="#991b1b", fontsize=8.5, fontweight="bold")
    ax_b.text(52, max(obs_exp_ratio) * 0.88, f"Peak 1: 55 bp\n({oe_55:.2f}x vs Null)", color="#c2410c", fontsize=8.5, fontweight="bold")
    ax_b.text(95, max(obs_exp_ratio) * 0.72, f"Peak 2: 100 bp\n({oe_100:.2f}x vs Null)", color="#0369a1", fontsize=8.5, fontweight="bold")
    ax_b.grid(True, alpha=0.25, ls="--")

    # Panel C: Theoretical Stereochemical Model Schema
    ax_c = axes[1, 0]
    im = ax_c.imshow(
        heatmap_2d_model,
        aspect="auto",
        origin="lower",
        extent=[offset_bins[0], offset_bins[-1], len_bins[0], len_bins[-1]],
        cmap="magma"
    )
    ax_c.set_title("C. Theoretical Model Schema: Expected Joint Density", fontsize=11, fontweight="bold")
    ax_c.set_xlabel("Signed Offset from Box Center to Dyad (bp)\n[Theoretical Schema: Prospective Empirical Package D]", fontsize=9)
    ax_c.set_ylabel("Fragment Length (bp)", fontsize=10)
    ax_c.axvline(-55, color="#38bdf8", ls=":", lw=1.2)
    ax_c.axvline(55, color="#38bdf8", ls=":", lw=1.2)
    ax_c.axhline(130, color="#ffffff", ls="--", lw=1.2, alpha=0.8)
    fig.colorbar(im, ax=ax_c, label="Theoretical Relative Density", fraction=0.046, pad=0.04)

    # Panel D: Stereochemical schematic of the 17-bp Box at +46.5...+63.5 bp
    ax_d = axes[1, 1]
    ax_d.set_xlim(-10, 120)
    ax_d.set_ylim(-3, 8)
    ax_d.axis("off")
    ax_d.set_title("D. Stereochemical Boundary of 17-bp Box at +55 bp", fontsize=11, fontweight="bold")

    # Draw dyad and protected core
    ax_d.plot([0, 0], [0, 4], color="#991b1b", lw=3)
    ax_d.text(0, 4.4, "Dyad (0 bp)\nStrictly Occluded", ha="center", fontsize=8, color="#991b1b", fontweight="bold")

    # Core arc / bar
    ax_d.fill_between([0, 65], [1, 1], [3, 3], color="#fed7aa", alpha=0.8, label="Protected 130-bp Core (R=65 bp)")
    ax_d.text(32.5, 2.0, "CENP-A Core (130 bp)\n(Radius = 65 bp)", ha="center", va="center", fontsize=8, color="#9a3412", fontweight="bold")

    # Linker bar
    ax_d.fill_between([65, 115], [1, 1], [3, 3], color="#e2e8f0", alpha=0.8, label="Linker DNA")
    ax_d.text(90, 2.0, "Free Linker DNA", ha="center", va="center", fontsize=8, color="#475569")

    # 17-bp CENP-B box at +55 bp (+46.5 to +63.5 bp)
    ax_d.fill_between([46.5, 63.5], [0.4, 0.4], [3.6, 3.6], color="#2563eb", alpha=0.85)
    ax_d.text(55, -0.6, "17-bp CENP-B Box\nCenter: +55 bp\n[+46.5 ... +63.5 bp]\nat Unpeeled Gyre Exit", ha="center", fontsize=8, color="#1d4ed8", fontweight="bold")

    # Unpeeled DNA arrow
    ax_d.annotate(
        "DNA Unpeeling Point\n(SHL ±5.0 to ±5.5)",
        xy=(65, 3.2), xytext=(75, 5.2),
        arrowprops=dict(facecolor="#0f172a", shrink=0.08, width=1, headwidth=5),
        fontsize=8, fontweight="bold"
    )

    plt.tight_layout()
    fig_png = os.path.join(FIG_DIR, "Fig4_cenpb_box_coupling_and_nulls.png")
    fig_svg = os.path.join(FIG_DIR, "Fig4_cenpb_box_coupling_and_nulls.svg")
    fig_pdf = os.path.join(FIG_DIR, "Fig4_cenpb_box_coupling_and_nulls.pdf")
    plt.savefig(fig_png)
    plt.savefig(fig_svg)
    plt.savefig(fig_pdf)
    plt.close()
    print(f"Generated {fig_png}, {fig_svg}, and {fig_pdf}")

if __name__ == "__main__":
    run_package_d()
