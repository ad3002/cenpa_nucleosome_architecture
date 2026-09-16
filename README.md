# Human CENP-A Nucleosome Core Architecture & Linker Geometry

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.placeholder.svg)](https://doi.org/10.5281/zenodo.placeholder)
[![Reproducibility](https://img.shields.io/badge/reproducibility-verified-brightgreen.svg)](run_reproduction.sh)
[![Assembly](https://img.shields.io/badge/genome-T2T--CHM13v2.0-blue.svg)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_009914755.4/)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)

This repository provides the complete, self-contained codebase, processed data tables, figures, and publication manuscript for our study:

> **"Human CENP-A Nucleosomes Form an Open 125–130 bp Particle Phased to a 340 bp Alpha-Satellite Dimer Lattice with Bipartite CENP-B Linker Geometry"**  
> *Aleksey Komissarov, Marina Popova, and Collaborators (2026).*

---

## Key Scientific Findings

1. **Native CENP-A Core Particle is 125–130 bp:**  
   Direct paired-end MNase ChIP-seq on human T2T centromeres (`SRR13278683`, 4.29M primary proper-pair fragments across 744 alpha-satellite arrays) reveals a dominant mononucleosome protection mode at **125–130 bp**. The ~80 bp hemisome (<2.1%) and the canonical 147 bp closed octamer (<0.3%) are quantitatively rejected *in vivo*, confirming that the terminal ~10 bp of DNA on each flank unpeel from the CENP-A octamer.
2. **Bipartite CENP-B Box Linker Positioning:**  
   Distance from the CENP-A dyad to 126,969 canonical 17-bp CENP-B boxes reveals severe dyad occlusion (0–15 bp, **>62-fold depletion**) and two sharp peaks:
   - **Peak 1 at 50–55 bp** (SHL $\pm 5.5$, directly at the unpeeled gyre exit boundary);
   - **Peak 2 at 85–100 bp** (free inter-nucleosomal linker DNA).
3. **Strict 340-bp Dimer Lattice Phasing in the CDR:**  
   Spatial autocorrelation (phasogram) of 411,919 dyads inside the hypomethylated Centromere Dip Region (CDR) demonstrates an alternating 150/190 bp monomer spacing (mean **170 bp**) that locks into a massive global maximum at **340 bp ($2 \times 170$ bp)** dinucleosome periodicity ($N = 761,698$ pairs).
4. **Bimodal Phase Transition via H1 Exclusion:**  
   In the heterochromatic periphery (Non-CDR), linker histone H1 and dense CpG methylation (85%) compress the repeat to **160 bp** (13 bp linkers). Inside the CDR, hypomethylation (25% 5mC) and unpeeled CENP-A ends evict H1, expanding linkers to 30–50 bp and generating an accessible 340-bp lattice for CENP-B dimer clamping.

---

## Repository Structure

```
.
├── README.md                      # This documentation & reproduction guide
├── run_reproduction.sh            # Master one-click reproduction script
├── paper/
│   ├── manuscript.md              # Full, watertight manuscript text
│   └── figures/
│       ├── Fig1_cenpa_core_and_box_geometry.png (.pdf, .svg)
│       └── Fig2_cdr_phasogram_and_chromatin_state.png (.pdf, .svg)
├── data/
│   ├── chm13_cdr_intervals.bed            # 23 CHM13 CDR genomic coordinates
│   ├── input_mnase_fragment_length_hist.tsv # Input MNase fragment distribution
│   ├── cenpa_chip_fragment_length_hist.tsv  # CENP-A ChIP fragment distribution
│   ├── input_box_to_dyad_distance.tsv     # Input dyad-to-box distance
│   ├── cenpa_box_to_dyad_distance.tsv     # CENP-A dyad-to-box distance
│   ├── cenpa_cdr_phasogram.tsv            # CENP-A CDR dyad autocorrelation
│   ├── cenpa_per_chromosome_summary.tsv   # Per-chromosome metrics across 23 chr
│   └── input_per_chromosome_summary.tsv   # Input per-chromosome metrics
└── scripts/
    ├── 01_fetch_and_align.sh      # Pipeline to fetch raw SRA slices and run BWA-MEM
    ├── 02_analyze_particles.py    # Extracts particle lengths, box distances from BAM
    ├── 03_compute_phasogram.py    # Spatial autocorrelation calculation
    ├── 04_plot_figures.py         # Generates vector publication figures (PDF/SVG/PNG)
    └── sync_paired.py             # FASTQ paired-end record synchronization
```

---

## Reproduction Guide

### Option A: Instant Quick Reproduction (< 30 seconds)
To regenerate all publication figures from the curated summary data tables:

```bash
# Clone the repository (once remote is added)
git clone <remote_url> cenpa_nucleosome_architecture
cd cenpa_nucleosome_architecture

# Run the master reproduction script in quick mode
./run_reproduction.sh
```
This will immediately generate `Fig1_cenpa_core_and_box_geometry.{png,pdf,svg}` and `Fig2_cdr_phasogram_and_chromatin_state.{png,pdf,svg}` in `paper/figures/`.

---

### Option B: Full Raw End-to-End Alignment (~1–2 hours)
To download raw FASTQ slices from the NCBI/EBI SRA (`SRR13278681` and `SRR13278683`), align them with BWA-MEM against CHM13 alpha arrays, and recompute all tables from scratch:

```bash
# Requirements: bwa, samtools, python3 (numpy, matplotlib)
./run_reproduction.sh --full-raw
```

---

## System Requirements

- **Operating System:** Linux / macOS (tested on Ubuntu 22.04 LTS and macOS Darwin)
- **Software Dependencies:**
  - Python $\ge 3.9$
  - `numpy`, `matplotlib`, `scipy`
  - `samtools` $\ge 1.15$ (for BAM parsing)
  - `bwa` $\ge 0.7.17$ (only for full-raw alignment mode)
  - `curl` (for raw FASTQ download)

---

## Citation & Contact

If you use the data, scripts, or structural models in this repository, please cite:
```bibtex
@article{komissarov2026cenpa,
  title={Human CENP-A Nucleosomes Form an Open 125--130 bp Particle Phased to a 340 bp Alpha-Satellite Dimer Lattice with Bipartite CENP-B Linker Geometry},
  author={Komissarov, Aleksey and Popova, Marina and others},
  journal={bioRxiv / In Review},
  year={2026}
}
```
For questions, contact Aleksey Komissarov (`akomissarov@...`).
