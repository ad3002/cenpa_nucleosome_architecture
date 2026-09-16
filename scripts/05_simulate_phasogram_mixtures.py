#!/usr/bin/env python3
"""
05_simulate_phasogram_mixtures.py

Demonstrates the mathematical counterexample to single-molecule lattice inference
from bulk spatial autocorrelation (phasograms).

Compares:
1. Model A (Alternating): Single molecules have strictly alternating 150 bp and 190 bp steps.
2. Model B (Register Mixture): Two cell populations each have uniform 340 bp repeats,
   shifted relative to each other by 150 bp.
   Neither population has alternating 150/190 bp steps, yet bulk autocorrelation
   yields the identical peak spectrum (150, 190, 340, 490, 530, 680 bp).
"""

import os
import csv
import collections
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def compute_autocorrelation(coords_list, max_lag=800, bin_size=5):
    counts = collections.Counter()
    for coords in coords_list:
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
    # 1. Model A: 50 independent fibers with alternating 150 / 190 steps
    fibers_alternating = []
    for f in range(50):
        pos = f * 17  # arbitrary random offset
        coords = [pos]
        for step in range(25):
            delta = 150 if step % 2 == 0 else 190
            pos += delta
            coords.append(pos)
        fibers_alternating.append(coords)

    # 2. Model B: Two cell populations with pure 340 bp spacing
    # Pop 1: 340 * k
    # Pop 2: 340 * k + 150
    # In neither population is there any 150/190 alternation on an individual fiber!
    fibers_pop1 = []
    for f in range(25):
        offset = f * 13
        coords = [offset + 340 * k for k in range(25)]
        fibers_pop1.append(coords)

    fibers_pop2 = []
    for f in range(25):
        offset = f * 13 + 150
        coords = [offset + 340 * k for k in range(25)]
        fibers_pop2.append(coords)

    fibers_mixture = fibers_pop1 + fibers_pop2

    # Bulk autocorrelation: in bulk sequencing, fragments from all fibers are pooled
    # pooled coordinates within a local array
    pooled_alternating = [c for fiber in fibers_alternating for c in fiber]
    pooled_mixture = [c for fiber in fibers_mixture for c in fiber]

    # Compute fiber-resolved vs bulk pooled autocorrelation
    # When bulk reads map to the same reference array, dyads are pooled into a single coordinate list
    bulk_alt_counts = compute_autocorrelation([pooled_alternating], max_lag=800, bin_size=5)
    bulk_mix_counts = compute_autocorrelation([pooled_mixture], max_lag=800, bin_size=5)

    # Export TSV
    out_dir = os.path.join(os.path.dirname(__file__), "..", "data")
    tsv_path = os.path.join(out_dir, "phasogram_simulation_comparison.tsv")
    with open(tsv_path, "w") as f:
        f.write("lag_bp\tmodel_a_alternating_bulk\tmodel_b_mixture_bulk\n")
        for lag in range(0, 805, 5):
            f.write(f"{lag}\t{bulk_alt_counts[lag]}\t{bulk_mix_counts[lag]}\n")
    print(f"Wrote {tsv_path}")

    # Plot figure
    fig_dir = os.path.join(os.path.dirname(__file__), "..", "paper", "figures")
    os.makedirs(fig_dir, exist_ok=True)
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True, dpi=300)

    lags = list(range(0, 805, 5))
    alt_vals = [bulk_alt_counts[l] for l in lags]
    mix_vals = [bulk_mix_counts[l] for l in lags]

    ax1.plot(lags, alt_vals, color="#c2410c", lw=2, label="Model A: Alternating Lattice (150-190-150-190)")
    ax1.set_title("Model A: Strict Intramolecular Alternation (150 / 190 bp steps)", fontsize=11, fontweight="bold")
    ax1.set_ylabel("Pairwise Autocorrelation", fontsize=10)
    ax1.grid(True, alpha=0.3, ls="--")
    for peak in [150, 190, 340, 490, 530, 680]:
        ax1.axvline(peak, color="#9a3412", ls=":", alpha=0.5)

    ax2.plot(lags, mix_vals, color="#0369a1", lw=2, label="Model B: Mixture of 2 Shifted Registers (both 340 bp)")
    ax2.set_title("Model B: Superposition of 2 Shifted 340-bp Registers (Phase Shift = 150 bp)", fontsize=11, fontweight="bold")
    ax2.set_xlabel("Inter-Dyad Distance (bp)", fontsize=10)
    ax2.set_ylabel("Pairwise Autocorrelation", fontsize=10)
    ax2.grid(True, alpha=0.3, ls="--")
    for peak in [150, 190, 340, 490, 530, 680]:
        ax2.axvline(peak, color="#0284c7", ls=":", alpha=0.5)

    plt.tight_layout()
    fig_png = os.path.join(fig_dir, "Fig3_phasogram_mixture_models.png")
    fig_svg = os.path.join(fig_dir, "Fig3_phasogram_mixture_models.svg")
    plt.savefig(fig_png)
    plt.savefig(fig_svg)
    plt.close()
    print(f"Generated {fig_png} and {fig_svg}")

if __name__ == "__main__":
    run_simulation()
