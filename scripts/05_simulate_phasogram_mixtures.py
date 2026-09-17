#!/usr/bin/env python3
"""
05_simulate_phasogram_mixtures.py

Demonstrates the mathematical equivalence between:
1. Model A (Alternating Lattice): An individual chromatin fiber with consecutive
   steps alternating between 150 bp and 190 bp (dyads at 0, 150, 340, 490, 680, 830, ...).
2. Model B (Register Mixture): Two cell subpopulations each possessing uniform 340-bp repeats,
   with Population 2 shifted by 150 bp relative to Population 1:
   - Pop 1: 0, 340, 680, 1020, ...
   - Pop 2: 150, 490, 830, 1170, ...
   When pooled in bulk sequencing at a single locus, the combined coordinate set
   is identical: {0, 150, 340, 490, 680, 830, ...}.

This proves algebraically that bulk spatial autocorrelation cannot distinguish
between true intramolecular alternation and a 50:50 mixture of two uniform 340-bp registers.
"""

import os
import csv
from collections import Counter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def compute_autocorrelation(coords, max_lag=800, bin_size=5):
    counts = Counter()
    sc = sorted(coords)
    n = len(sc)
    for i in range(n):
        for j in range(i + 1, n):
            lag = sc[j] - sc[i]
            if lag > max_lag:
                break
            b = int(round(lag / bin_size) * bin_size)
            counts[b] += 1
    return counts

def run_simulation():
    K = 30  # number of repeat units
    
    # Model A: Single fiber with alternating 150 bp and 190 bp steps
    coords_a = [0]
    pos = 0
    for step in range(K * 2):
        delta = 150 if step % 2 == 0 else 190
        pos += delta
        coords_a.append(pos)
        
    # Model B: Two independent cell populations with uniform 340 bp spacing
    pop1 = [340 * k for k in range(K + 1)]
    pop2 = [340 * k + 150 for k in range(K)]
    coords_b = sorted(pop1 + pop2)

    # In bulk sequencing, reads mapping to this array sample both populations.
    # Therefore, the pooled set of observed dyad positions across the population is coords_b.
    counts_a = compute_autocorrelation(coords_a, max_lag=800, bin_size=5)
    counts_b = compute_autocorrelation(coords_b, max_lag=800, bin_size=5)

    lags = list(range(0, 805, 5))
    vals_a = [counts_a[l] for l in lags]
    vals_b = [counts_b[l] for l in lags]

    # Verification: vals_a and vals_b must be IDENTICAL because coords_a == coords_b
    assert coords_a == coords_b, "Mathematical coordinates must be identical!"
    assert vals_a == vals_b, "Autocorrelation spectrum must be identically equal!"

    # Save TSV
    out_dir = os.path.join(os.path.dirname(__file__), "..", "data")
    tsv_path = os.path.join(out_dir, "phasogram_simulation_comparison.tsv")
    with open(tsv_path, "w") as f:
        f.write("lag_bp\tmodel_a_alternating\tmodel_b_register_mixture\tdifference\n")
        for l, a, b in zip(lags, vals_a, vals_b):
            f.write(f"{l}\t{a}\t{b}\t{a - b}\n")
    print(f"Wrote {tsv_path} (Identical proof verified: max difference = 0)")

    # Plot Figure 3
    fig_dir = os.path.join(os.path.dirname(__file__), "..", "paper", "figures")
    os.makedirs(fig_dir, exist_ok=True)
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(9, 7.5), sharex=True, dpi=300)

    # Panel 1: Model A
    ax1.plot(lags, vals_a, color="#c2410c", lw=2, label="Model A: Intramolecular Alternation (150 / 190 bp)")
    ax1.set_title("Model A: Single Molecule Alternating Lattice (Steps: 150 bp, 190 bp, 150 bp, ...)", fontsize=10, fontweight="bold")
    ax1.set_ylabel("Pair Count", fontsize=9)
    ax1.grid(True, alpha=0.25, ls="--")
    for peak in [150, 190, 340, 490, 530, 680]:
        ax1.axvline(peak, color="#9a3412", ls=":", alpha=0.6)
        ax1.text(peak, max(vals_a)*0.88, f"{peak}", color="#9a3412", fontsize=7.5, ha="center", rotation=90)
    ax1.legend(loc="upper right", fontsize=8)

    # Panel 2: Model B
    ax2.plot(lags, vals_b, color="#0369a1", lw=2, label="Model B: Mixture of 2 Shifted 340-bp Registers (Phase Shift = 150 bp)")
    ax2.set_title("Model B: Ensemble Mixture of 2 Uniform 340-bp Registers (Pop 1: 340k; Pop 2: 340k + 150)", fontsize=10, fontweight="bold")
    ax2.set_ylabel("Pair Count", fontsize=9)
    ax2.grid(True, alpha=0.25, ls="--")
    for peak in [150, 190, 340, 490, 530, 680]:
        ax2.axvline(peak, color="#0284c7", ls=":", alpha=0.6)
        ax2.text(peak, max(vals_b)*0.88, f"{peak}", color="#0284c7", fontsize=7.5, ha="center", rotation=90)
    ax2.legend(loc="upper right", fontsize=8)

    # Panel 3: Direct Overlay & Residual Difference
    ax3.plot(lags, [a - b for a, b in zip(vals_a, vals_b)], color="#15803d", lw=1.5, label="Residual Difference (Model A − Model B)")
    ax3.set_title("Mathematical Equivalence: Residual Difference is Exactly Zero (r = 1.000)", fontsize=10, fontweight="bold")
    ax3.set_xlabel("Inter-Dyad Distance (bp)", fontsize=10)
    ax3.set_ylabel("Residual (A − B)", fontsize=9)
    ax3.set_ylim(-5, 5)
    ax3.grid(True, alpha=0.25, ls="--")
    ax3.legend(loc="upper right", fontsize=8)

    plt.tight_layout()
    fig_png = os.path.join(fig_dir, "Fig3_phasogram_mixture_models.png")
    fig_svg = os.path.join(fig_dir, "Fig3_phasogram_mixture_models.svg")
    fig_pdf = os.path.join(fig_dir, "Fig3_phasogram_mixture_models.pdf")
    plt.savefig(fig_png)
    plt.savefig(fig_svg)
    plt.savefig(fig_pdf)
    plt.close()
    print(f"Generated {fig_png}, {fig_svg}, and {fig_pdf}")

if __name__ == "__main__":
    run_simulation()
