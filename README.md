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

1. **Native CENP-A Open Core Particle (125–133 bp):**  
   Direct paired-end MNase ChIP-seq on human T2T centromeres (`SRR13278683`, 4.29M primary proper-pair fragments across 744 alpha-satellite arrays) reveals a single-base mononucleosome protection mode at **133 bp** ($N = 212,205$ fragments globally), with $177,473$ fragments at 130 bp and **84.29%** of all fragments ($N = 3,616,490$) concentrated in the 110–140 bp window (with 130 bp representing the modal 5-bp bin). In this library, the canonical 150-bp octamer is depleted **177.28-fold** relative to the 133-bp mode ($N = 1,197$, 0.0279%; 148.26-fold relative to 130 bp), and sub-85 bp fragments represent 1.53%, confirming that terminal gyres unpeel from the CENP-A octamer core in native centromeric chromatin.
2. **Bipartite CENP-B Box Architecture:**  
   Distance from the CENP-A dyad to 126,969 canonical 17-bp CENP-B boxes reveals severe dyad occlusion (0–15 bp, **66.0-fold contrast** relative to Peak 1; 15.36-fold depletion relative to an empirical stepwise null model baseline) and two major peaks:
   - **Peak 1 at 55-bp bin** (SHL $\pm 5.0\text{--}5.5$, $[55, 60)$ bp with box centered at +55 bp; 4.30x enrichment vs stepwise null baseline);
   - **Peak 2 at 85–100 bp** (free inter-nucleosomal linker DNA; 6.54x enrichment vs stepwise null baseline at 100 bp).
3. **Dominant 340-bp Dimer Periodicity in the CDR:**  
   Spatial autocorrelation (phasogram) of 411,919 dyads inside the hypomethylated Centromere Dip Region (CDR) exhibits bimodal monomer lags at 150 bp and 190 bp (mean **170 bp**) that lock into a dominant non-zero maximum at **340 bp ($2 \times 170$ bp)** dinucleosome periodicity ($N = 761,698$ pairs in CDR; 591,711 in Non-CDR).
4. **Ensemble Phasing & Chromatin Architecture Models:**  
   Simulating bulk autocorrelation proves that an intramolecular alternating lattice (150/190 bp steps) and a 50:50 mixture of two uniform 340-bp registers produce identically matching autocorrelation spectra ($r = 1.000$, residual difference = 0). Expanded 20 bp and 60 bp linkers (mean 40 bp) and linker histone H1 exclusion represent mechanistic models inferred from bulk spacing and terminal unpeeling, providing an accessible chromatin framework for CENP-B dimer cross-linking.

---

## Repository Structure

```
.
├── README.md                          # Documentation & reproduction guide
├── run_reproduction.sh                # Master verification suite with dynamic assert checks (<4s)
├── paper/
│   ├── manuscript.md                  # Comprehensive manuscript text
│   ├── index.html                     # Interactive publication article
│   └── figures/
│       ├── Fig1_cenpa_core_and_box_geometry.{png,pdf,svg}
│       ├── Fig2_cdr_phasogram_and_chromatin_state.{png,pdf,svg}
│       ├── Fig3_phasogram_mixture_models.{png,pdf,svg}
│       └── Fig4_cenpb_box_coupling_and_nulls.{png,pdf,svg}
├── data/
│   ├── ledger_manifest.tsv            # Dynamically generated single-source-of-truth ledger
│   ├── metrics.json                   # Machine-readable verified parameters
│   ├── replicate_metadata_manifest.tsv# Authenticated cross-cohort replicate manifest
│   ├── chm13_cdr_intervals.bed        # 23 CHM13 CDR genomic coordinates
│   ├── input_mnase_fragment_length_hist.tsv # Input MNase fragment distribution
│   ├── cenpa_chip_fragment_length_hist.tsv  # CENP-A ChIP fragment distribution
│   ├── cenpa_box_to_dyad_distance.tsv # CENP-A dyad-to-box distance
│   ├── cenpa_box_directional_and_nulls.tsv  # Dyad-to-box vs Stepwise Null calibration
│   ├── cenpa_cdr_phasogram.tsv        # CENP-A CDR dyad autocorrelation
│   ├── phasogram_simulation_comparison.tsv  # Algebraic equivalence of mixture models
│   ├── cenpa_per_chromosome_summary.tsv     # Per-chromosome metrics across 23 chr
│   └── archive/                       # Archived prospective/synthetic scratch tables
└── scripts/
    ├── generate_ledger.py             # Dynamically reconciles ledger sums, modes, and ratios
    ├── build_replicate_manifest.py    # Generates authenticated public SRA metadata
    ├── sync_paired.py                 # Streaming paired-end FASTQ synchronizer & validator
    ├── annotate_cenpb_boxes.py        # Exact 17-bp CENP-B box motif annotator across alpha arrays
    ├── 01_fetch_and_align.sh          # SRA fetch and BWA-MEM alignment
    ├── 02_analyze_particles.py        # Particle lengths & box distances from BAM
    ├── 03_compute_phasogram.py        # Distance-bounded spatial autocorrelation
    ├── 04_plot_figures.py             # Primary empirical figures (Fig 1 & 2)
    ├── 05_simulate_phasogram_mixtures.py    # Register mixture simulation (Fig 3)
    └── 06_analyze_box_coupling_and_nulls.py # CENP-B coupling & Stepwise Null Baseline (Fig 4)
```

---

## Reproduction Guide

### Option A: Instant Quick Reproduction (< 4 seconds)
To verify all dynamic ledger assertions and regenerate all 4 publication figures from the authentic data tables:

```bash
# Clone the repository
git clone https://github.com/ad3002/cenpa_nucleosome_architecture.git
cd cenpa_nucleosome_architecture

# Run the master reproduction script
bash run_reproduction.sh
```
This runs independent dynamic recalculations from raw tables and regenerates:
- `Fig1_cenpa_core_and_box_geometry.{png,pdf,svg}`
- `Fig2_cdr_phasogram_and_chromatin_state.{png,pdf,svg}`
- `Fig3_phasogram_mixture_models.{png,pdf,svg}`
- `Fig4_cenpb_box_coupling_and_nulls.{png,pdf,svg}`
in `paper/figures/`.

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
