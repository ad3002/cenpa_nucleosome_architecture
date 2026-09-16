#!/usr/bin/env python3
"""
04_plot_figures.py

Generates publication-quality figures for the CENP-A nucleosome architecture paper.
Outputs:
- paper/figures/Fig1_cenpa_core_and_box_geometry.png (.pdf, .svg)
- paper/figures/Fig2_cdr_phasogram_and_chromatin_state.png (.pdf, .svg)
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "paper", "figures")
os.makedirs(OUT_DIR, exist_ok=True)

# Set global publication styling
plt.rcParams["font.sans-serif"] = "Helvetica", "Arial", "DejaVu Sans"
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["figure.dpi"] = 300
plt.rcParams["axes.linewidth"] = 1.0
plt.rcParams["xtick.major.width"] = 1.0
plt.rcParams["ytick.major.width"] = 1.0
plt.rcParams["xtick.minor.width"] = 0.6
plt.rcParams["ytick.minor.width"] = 0.6
plt.rcParams["axes.titlesize"] = 11
plt.rcParams["axes.labelsize"] = 10
plt.rcParams["xtick.labelsize"] = 9
plt.rcParams["ytick.labelsize"] = 9
plt.rcParams["legend.fontsize"] = 8.5

# ==============================================================================
# FIGURE 1: CENP-A Open Core Particle and Linker Architecture
# ==============================================================================
fig1 = plt.figure(figsize=(10.5, 4.2))
gs1 = gridspec.GridSpec(1, 2, width_ratios=[1.1, 1.1], wspace=0.28)

# Panel A: Fragment length distribution (Input vs CENP-A ChIP)
ax1 = fig1.add_subplot(gs1[0])

input_hist = os.path.join(DATA_DIR, "input_mnase_fragment_length_hist.tsv")
cenpa_hist = os.path.join(DATA_DIR, "cenpa_chip_fragment_length_hist.tsv")

inp_data = np.loadtxt(input_hist, skiprows=1)
cen_data = np.loadtxt(cenpa_hist, skiprows=1)

# Normalize distributions to percentage in 60-200 bp range
mask_inp = (inp_data[:, 0] >= 60) & (inp_data[:, 0] <= 200)
mask_cen = (cen_data[:, 0] >= 60) & (cen_data[:, 0] <= 200)

x_inp = inp_data[mask_inp, 0]
y_inp = inp_data[mask_inp, 1]
y_inp_norm = 100.0 * y_inp / np.sum(y_inp)

x_cen = cen_data[mask_cen, 0]
y_cen = cen_data[mask_cen, 1]
y_cen_norm = 100.0 * y_cen / np.sum(y_cen)

ax1.plot(x_inp, y_inp_norm, color="#555555", lw=2.2, label="Input MNase (Total Centromeric H3)", alpha=0.9)
ax1.fill_between(x_inp, y_inp_norm, color="#888888", alpha=0.15)

ax1.plot(x_cen, y_cen_norm, color="#D9381E", lw=2.5, label="CENP-A MNase ChIP-seq (Native kinetochore)")
ax1.fill_between(x_cen, y_cen_norm, color="#D9381E", alpha=0.20)

# Annotations
ax1.axvline(147, color="#333333", linestyle="--", lw=1.0, alpha=0.7)
ax1.text(149, np.max(y_inp_norm) * 0.92, "Canonical H3 core\n(147 bp)", color="#333333", fontsize=8, ha="left")

ax1.axvline(130, color="#D9381E", linestyle="--", lw=1.2, alpha=0.9)
ax1.text(128, np.max(y_cen_norm) * 0.95, "CENP-A open core\n(125–130 bp)", color="#D9381E", fontsize=8.5, fontweight="bold", ha="right")

ax1.axvline(80, color="#2B6CB0", linestyle=":", lw=1.0, alpha=0.6)
ax1.text(82, 1.0, "Hemisome null\n(~80 bp: <2%)", color="#2B6CB0", fontsize=7.5, ha="left")

ax1.set_xlabel("Fragment Length (bp)")
ax1.set_ylabel("Normalized Fragment Density (%)")
ax1.set_title("A. Native CENP-A Nucleosome Core Particle Protection", loc="left", fontweight="bold")
ax1.set_xlim(60, 200)
ax1.set_ylim(0, max(np.max(y_inp_norm), np.max(y_cen_norm)) * 1.12)
ax1.legend(loc="upper left", frameon=True, framealpha=0.9)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

# Panel B: Dyad to CENP-B box distance distribution
ax2 = fig1.add_subplot(gs1[1])

cenpa_dist_file = os.path.join(DATA_DIR, "cenpa_box_to_dyad_distance.tsv")
dist_data = np.loadtxt(cenpa_dist_file, skiprows=1)

d_x = dist_data[:, 0]
d_y = dist_data[:, 1] / 1000.0  # in thousands

# Mask up to 200 bp
m_dist = d_x <= 180
d_x = d_x[m_dist]
d_y = d_y[m_dist]

ax2.plot(d_x, d_y, color="#2B6CB0", lw=2.4, marker="o", markersize=3.5, label="CENP-A Dyad to CENP-B Box")
ax2.fill_between(d_x, d_y, color="#2B6CB0", alpha=0.18)

# Core particle boundary
ax2.axvspan(0, 65, color="#D9381E", alpha=0.08, label="CENP-A Core Radius (65 bp)")
ax2.axvspan(65, 180, color="#2B6CB0", alpha=0.05, label="Inter-nucleosomal Linker DNA")

# Peak labels
ax2.annotate("Peak 1: 55 bp\n(Gyre Exit / SHL ±5.5)\n164.7k events", 
             xy=(55, 164.7), xytext=(20, 135),
             arrowprops=dict(facecolor="#D9381E", shrink=0.08, width=1, headwidth=5),
             fontsize=8, fontweight="bold", color="#D9381E")

ax2.annotate("Peak 2: 90–100 bp\n(Free Linker DNA)\n125.4k events", 
             xy=(95, 125.4), xytext=(115, 120),
             arrowprops=dict(facecolor="#2B6CB0", shrink=0.08, width=1, headwidth=5),
             fontsize=8, fontweight="bold", color="#2B6CB0")

# Dyad axis depletion
ax2.annotate("Dyad Occlusion\n(0–15 bp: >62x depletion)", 
             xy=(5, 5), xytext=(10, 45),
             arrowprops=dict(facecolor="#333333", shrink=0.08, width=1, headwidth=4),
             fontsize=7.5, color="#444444")

ax2.set_xlabel("Distance from CENP-A Dyad to CENP-B Box (bp)")
ax2.set_ylabel("Fragment Midpoint Count (x 1,000)")
ax2.set_title("B. Bipartite CENP-B Box Linker Positioning", loc="left", fontweight="bold")
ax2.set_xlim(0, 180)
ax2.set_ylim(0, 185)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.legend(loc="upper right", frameon=True, framealpha=0.9)

for ext in ["png", "pdf", "svg"]:
    fig1.savefig(os.path.join(OUT_DIR, f"Fig1_cenpa_core_and_box_geometry.{ext}"), dpi=300, bbox_inches="tight")
plt.close(fig1)

print("Saved Figure 1 (PNG, PDF, SVG)")

# ==============================================================================
# FIGURE 2: Centromeric Dip Region Phasogram & Dimer Lattice
# ==============================================================================
fig2 = plt.figure(figsize=(10.5, 4.2))
gs2 = gridspec.GridSpec(1, 2, width_ratios=[1.15, 0.95], wspace=0.28)

# Panel A: CENP-A CDR Phasogram (Autocorrelation)
ax3 = fig2.add_subplot(gs2[0])

phas_file = os.path.join(DATA_DIR, "cenpa_cdr_phasogram.tsv")
phas_data = np.loadtxt(phas_file, skiprows=1)

p_dist = phas_data[:, 0]
p_cdr = phas_data[:, 1] / 1000.0      # thousands in CDR
p_noncdr = phas_data[:, 2] / 1000.0   # thousands in Non-CDR

m_p = (p_dist >= 100) & (p_dist <= 500)
p_dist_sub = p_dist[m_p]
p_cdr_sub = p_cdr[m_p]
p_noncdr_sub = p_noncdr[m_p]

ax3.plot(p_dist_sub, p_cdr_sub, color="#D9381E", lw=2.2, label="CENP-A Dyads inside CDR (Active Kinetochore)")
ax3.fill_between(p_dist_sub, p_cdr_sub, color="#D9381E", alpha=0.15)

# Periodic lattice markers
ax3.axvline(170, color="#666666", linestyle=":", lw=1.0)
ax3.text(172, 540, "171 bp\n(1x Alpha monomer)", color="#444444", fontsize=7.5)

ax3.axvline(340, color="#D9381E", linestyle="--", lw=1.3)
ax3.annotate("340 bp Dimer Lattice\n(2 x 170 bp: 761.7k pairs)", 
             xy=(340, 761.7), xytext=(360, 680),
             arrowprops=dict(facecolor="#D9381E", shrink=0.08, width=1, headwidth=5),
             fontsize=8.5, fontweight="bold", color="#D9381E")

# Split peaks at monomer
ax3.plot(150, p_cdr_sub[p_dist_sub == 150], "o", color="#2B6CB0", markersize=6)
ax3.text(142, 550, "150 bp", color="#2B6CB0", fontsize=8, fontweight="bold", ha="right")

ax3.plot(190, p_cdr_sub[p_dist_sub == 190], "o", color="#2B6CB0", markersize=6)
ax3.text(198, 500, "190 bp", color="#2B6CB0", fontsize=8, fontweight="bold", ha="left")

ax3.set_xlabel("Distance Between Adjacent CENP-A Dyads (bp)")
ax3.set_ylabel("Pairwise Spatial Dyad Pairs (x 1,000)")
ax3.set_title("A. CENP-A Spatial Autocorrelation (CDR Phasogram)", loc="left", fontweight="bold")
ax3.set_xlim(100, 480)
ax3.set_ylim(0, 850)
ax3.spines["top"].set_visible(False)
ax3.spines["right"].set_visible(False)
ax3.legend(loc="upper left", frameon=True, framealpha=0.9)

# Panel B: Summary of Chromatin State Transition
ax4 = fig2.add_subplot(gs2[1])

categories = ["Mononucleosome Core", "Linker DNA Length", "Nucleosome Repeat (NRL)", "Linker Histone H1", "CENP-B Box Status"]
periphery_vals = [147, 13, 160, 1.0, 0.2]
cdr_vals = [128, 42, 170, 0.0, 1.0]

# Bar comparison of physical lengths
labels = ["Core (bp)", "Linker (bp)", "NRL (bp)"]
noncdr_lens = [147, 13, 160]
cdr_lens = [128, 42, 170]

x = np.arange(len(labels))
width = 0.35

rects1 = ax4.bar(x - width/2, noncdr_lens, width, label="Periphery (Non-CDR, H3K9me3+)", color="#555555", alpha=0.85)
rects2 = ax4.bar(x + width/2, cdr_lens, width, label="Kinetochore (CDR, CENP-A+)", color="#D9381E", alpha=0.85)

# Value annotations on bars
for rect in rects1:
    h = rect.get_height()
    ax4.text(rect.get_x() + rect.get_width()/2., h + 3, f"{int(h)}", ha="center", va="bottom", fontsize=8, color="#333333")

for rect in rects2:
    h = rect.get_height()
    ax4.text(rect.get_x() + rect.get_width()/2., h + 3, f"{int(h)}", ha="center", va="bottom", fontsize=8, fontweight="bold", color="#D9381E")

ax4.set_ylabel("DNA Length (Base Pairs)")
ax4.set_title("B. Bimodal Centromeric Phase Transition", loc="left", fontweight="bold")
ax4.set_xticks(x)
ax4.set_xticklabels(labels, fontweight="bold")
ax4.set_ylim(0, 210)
ax4.spines["top"].set_visible(False)
ax4.spines["right"].set_visible(False)
ax4.legend(loc="upper left", frameon=True, framealpha=0.9)

# Additional state box
text_box = (
    "Chromatin State Partitioning:\n"
    "• Periphery: 5mC (85%), H1 bound,\n"
    "  dense 160 bp lattice, closed linkers.\n"
    "• CDR Core: 5mC (25%), H1 excluded,\n"
    "  open 128 bp core, 340 bp dimer lattice."
)
ax4.text(0.5, 0.52, text_box, transform=ax4.transAxes, fontsize=7.5,
         verticalalignment="top", bbox=dict(boxstyle="round,pad=0.5", facecolor="#FFF9E6", edgecolor="#E2C974", alpha=0.9))

for ext in ["png", "pdf", "svg"]:
    fig2.savefig(os.path.join(OUT_DIR, f"Fig2_cdr_phasogram_and_chromatin_state.{ext}"), dpi=300, bbox_inches="tight")
plt.close(fig2)

print("Saved Figure 2 (PNG, PDF, SVG)")
