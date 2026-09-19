# Human CENP-A Nucleosome Core Architecture & Linker Geometry

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.placeholder.svg)](https://doi.org/10.5281/zenodo.placeholder)
[![Reproducibility](https://img.shields.io/badge/reproducibility-verified-brightgreen.svg)](run_reproduction.sh)
[![Assembly](https://img.shields.io/badge/genome-T2T--CHM13v2.0-blue.svg)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_009914755.4/)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)

This repository provides the complete, self-contained codebase, processed data tables, figures, publication manuscript, and standalone reviewer verification package for:

> **"Human CENP-A Nucleosomes Form an Open 125–130 bp Particle Phased to a 340 bp Alpha-Satellite Dimer Lattice with Bipartite CENP-B Linker Geometry"**  
> *Aleksey Komissarov, Marina Popova, and Collaborators (2026).*

---

## Reviewer & Auditor Inspection Package

For peer reviewers, referees, and editors, a dedicated standalone package is prepared:
* **Reviewer Guide & Audit Traceability Matrix:** [`REVIEWER_PACKAGE_GUIDE.md`](REVIEWER_PACKAGE_GUIDE.md)
* **Standalone Distribution Archives:**
  * [`cenpa_reviewer_package_v3.1.tar.gz`](cenpa_reviewer_package_v3.1.tar.gz) (7.2 MB)
  * [`cenpa_reviewer_package_v3.1.zip`](cenpa_reviewer_package_v3.1.zip) (7.2 MB)
* **Compiled Publication Manuscripts:**
  * **Integrated Manuscript with Embedded Figures (PDF, 15 pages):** [`paper/manuscript_integrated.pdf`](paper/manuscript_integrated.pdf)
  * **Standard Submission Manuscript (PDF, 13 pages):** [`paper/manuscript.pdf`](paper/manuscript.pdf)
  * **Interactive Web Article (KaTeX + Sidebar):** [`paper/index.html`](paper/index.html)
  * **Markdown Source:** [`paper/manuscript.md`](paper/manuscript.md)
* **Supplementary Data Workbook (Excel, 8 Formatted Tabs):**
  * [`data/Supplementary_Tables_S1_to_S8.xlsx`](data/Supplementary_Tables_S1_to_S8.xlsx)

---

## Key Scientific Findings (Packages A–G)

1. **Native CENP-A Open Core Particle (125–133 bp):**  
   Direct paired-end MNase ChIP-seq on human T2T centromeres (`SRR13278683`, 4,290,331 primary proper-pair fragments across 744 alpha-satellite arrays) reveals a single-base mononucleosome protection mode at **133 bp** ($N = 212,205$ fragments globally), with $177,473$ fragments at 130 bp and **84.29%** of all fragments ($N = 3,616,490$) concentrated in the 110–140 bp window (with 130 bp representing the modal 5-bp bin). In this library, the canonical 150-bp octamer is depleted **177.28-fold** relative to the 133-bp mode ($N = 1,197$, 0.0279%; 148.26-fold relative to 130 bp), and sub-85 bp fragments represent 1.53%, confirming that terminal gyres unpeel from the CENP-A octamer core in native centromeric chromatin (**Fig. 1**).

2. **Reference-Free Physical Caliper & MAPQ Invariance (Packages B & C):**  
   - **Physical Caliper:** Direct measurement of physical insert length from adapter-bounded raw FASTQ read overlaps without reference alignment confirms the **133-bp mode** (88.18% in [110, 140] bp within $L \le 138$ bp) with **98.49%** ($N = 74,763 / 75,911$) exact base concordance ($R^2 = 0.8832$, weighted mean difference = -0.31 bp, median difference = 0.0 bp) with BWA-MEM BAM alignment (**Fig. 5A–C**).
   - **MAPQ Invariance:** Both multi-mapping fragments ($\text{MAPQ} = 0$, $N = 80,957$) and uniquely placed fragments ($\text{MAPQ} \ge 20$, $N = 2,942$) exhibit identical single-base protection modes at **133 bp** ($\Delta = 0$ bp), demonstrating complete immunity to repeat mapping ambiguity (**Fig. 5D**).

3. **Bipartite CENP-B Box Architecture (Package D):**  
   Distance from the CENP-A dyad to 126,969 canonical 17-bp CENP-B boxes reveals severe dyad occlusion (0–15 bp, **66.0-fold contrast** relative to Peak 1; 15.36-fold depletion relative to an empirical stepwise null baseline) and two major functional zones:
   - **Peak 1 at 55-bp bin** (SHL $\pm 5.0\text{--}5.5$, unpeeled gyre boundary; 4.30x enrichment vs stepwise null baseline);
   - **Peak 2 at 85–100 bp** (free inter-nucleosomal linker DNA; 6.54x enrichment vs stepwise null baseline at 100 bp) (**Fig. 4**).

4. **Dominant 340-bp Dimer Periodicity & Mathematical Equivalence (Package E):**  
   Spatial autocorrelation (phasogram) of 411,919 dyads inside the hypomethylated Centromere Dip Region (CDR) exhibits bimodal monomer lags at 150 bp and 190 bp (mean **170 bp**) that lock into a dominant non-zero maximum at **340 bp ($2 \times 170$ bp)** dinucleosome periodicity ($N = 761,698$ pairs in CDR; 591,711 in Non-CDR) (**Fig. 2**). Simulating bulk autocorrelation proves that an intramolecular alternating lattice (150/190 bp steps) and a 50:50 mixture of two uniform 340-bp registers produce identically matching autocorrelation spectra ($r = 1.000$, residual difference $\equiv 0$), defining an epistemic boundary requiring single-molecule long-read sequencing (**Fig. 3**).

5. **Local Intra-Array Epigenetic Contrast (Package F):**  
   Paired contrast within the **exact same continuous HOR array** across chromosomes (`hor_1_5`, `hor_8_2`, `hor_11_3`, `hor_X_1`) demonstrates **3.84-fold pooled CENP-A enrichment** in active CDR cores (4.374 rp/kb) over intra-array flanks (1.140 rp/kb; exact two-sided sign-test $p = 2.38 \times 10^{-7}$, two-sided Wilcoxon $p = 2.70 \times 10^{-5}$; chr11: 4.23 vs 0.85 rp/kb [4.98x]). Flanking 160-bp repeats leave ~13-bp linkers that sterically clash with the 17-bp CENP-B box and bind linker histone H1, whereas CDR 340-bp units provide 20/60 bp linkers that accommodate CENP-B boxes and exclude H1 (**Fig. 6**).

6. **Cross-Lineage Biological Replication (Package G):**  
   Biological validation across independent cohorts, karyotypes, and platforms confirms architectural conservation (**Fig. 7**, **Table 2**):
   - **CHM13 Biological Replicate 1 (`SRR13278684`):** Mode at **133 bp** (76.64% in core gate; 0.215% 150 bp, 26.91-fold depleted), FASTQ caliper mode at **133 bp**, and independent replication of the **340-bp dimer lattice peak** ($N = 779$ pairs in Rep 1; $N = 844$ in Rep 2).
   - **HG002 T2T Diploid (`SRR15395857`):** Diploid B-lymphoblastoid cell line (GM24385) exhibits conditional mononucleosome mode at **120 bp** (modal 5-bp bin at 125 bp), unconditioned mode 20 bp, caliper mode 90 bp, and canonical 150-bp depletion to **0.435%**.
   - **RPE-1 Non-Transformed Diploid (`SRR9201843`):** CENP-A mode at **175 bp**; canonical 150-bp octamers depleted to **0.608%** (<1% across all cohorts).
   - **RPE-1 CENP-B Architectural Comparator (`SRR9201844`):** Broad multi-protein kinetochore factor footprint (modal bin at 165 bp; tied single-base modes 126, 162, 163, 165 bp; only 0.74% in 45–65 bp), distinguishing sequence-specific complexes from nucleosome wrapping.

---

## Repository Structure

```text
.
├── README.md                          # Repository documentation & guide
├── REVIEWER_PACKAGE_GUIDE.md          # Comprehensive Reviewer & Auditor Guide
├── cenpa_reviewer_package_v3.1.tar.gz # Compressed reviewer distribution archive
├── cenpa_reviewer_package_v3.1.zip    # Compressed reviewer distribution archive (ZIP)
├── run_reproduction.sh                # Master verification suite with dynamic assert checks (~15s)
├── paper/
│   ├── manuscript_integrated.pdf      # Publication PDF with embedded high-resolution figures
│   ├── manuscript.pdf                 # Standard submission PDF manuscript
│   ├── manuscript.md                  # Markdown manuscript source
│   ├── index.html                     # Interactive web article
│   └── figures/
│       ├── Fig1_cenpa_core_and_box_geometry.{png,pdf,svg}
│       ├── Fig2_cdr_phasogram_and_chromatin_state.{png,pdf,svg}
│       ├── Fig3_phasogram_mixture_models.{png,pdf,svg}
│       ├── Fig4_cenpb_box_coupling_and_nulls.{png,pdf,svg}
│       ├── Fig5_physical_caliper_and_mapq_invariance.{png,pdf,svg}
│       ├── Fig6_intra_array_epigenetic_contrast.{png,pdf,svg}
│       └── Fig7_cross_lineage_replication.{png,pdf,svg}
├── data/
│   ├── Supplementary_Tables_S1_to_S8.xlsx # Formatted multi-tab Excel workbook
│   ├── ledger_manifest.tsv            # Dynamically generated single-source-of-truth ledger
│   ├── metrics.json                   # Machine-readable verified parameters
│   ├── replicate_metadata_manifest.tsv# Authenticated cross-cohort replicate manifest
│   ├── chm13_cdr_intervals.bed        # 23 CHM13 CDR genomic coordinates
│   ├── cenpa_chip_fragment_length_hist.tsv # CENP-A ChIP fragment distribution
│   ├── input_mnase_fragment_length_hist.tsv# Input MNase fragment distribution
│   ├── cenpa_cdr_phasogram.tsv        # CENP-A CDR dyad autocorrelation
│   ├── cenpa_per_chromosome_summary.tsv    # Per-chromosome metrics across 23 chr
│   ├── cenpa_box_directional_and_nulls.tsv # Dyad-to-box vs Stepwise Null calibration
│   ├── read_overlap_caliper_hist.tsv  # FASTQ physical read overlap distribution
│   ├── fragment_length_by_mapq.tsv    # Fragment lengths stratified by MAPQ
│   ├── caliper_vs_tlen_concordance.tsv# Concordance between Caliper and BAM TLEN
│   ├── intra_array_cdr_vs_flank_metrics.tsv# Intra-array CDR vs Flank read metrics
│   ├── intra_array_transition_summary.tsv  # Summary of intra-array density transition
│   ├── cross_lineage_metrics_summary.tsv   # Cross-lineage replication summary
│   ├── cross_lineage_length_distributions.tsv # Cross-lineage fragment length distributions
│   └── cross_lineage_phasograms.tsv   # Cross-lineage dyad phasograms
└── scripts/
    ├── build_reviewer_package.py      # Automated compiler for reviewer bundle & Excel tables
    ├── generate_ledger.py             # Dynamically reconciles ledger sums, modes, and ratios
    ├── build_replicate_manifest.py    # Generates authenticated public SRA metadata
    ├── sync_paired.py                 # Streaming paired-end FASTQ synchronizer & validator
    ├── annotate_cenpb_boxes.py        # Exact 17-bp CENP-B box motif annotator
    ├── 01_fetch_and_align.sh          # SRA fetch and BWA-MEM alignment
    ├── 02_analyze_particles.py        # Particle lengths & box distances from BAM
    ├── 03_compute_phasogram.py        # Distance-bounded spatial autocorrelation
    ├── 04_plot_figures.py             # Primary empirical figures (Fig 1 & 2)
    ├── 05_simulate_phasogram_mixtures.py    # Register mixture simulation (Fig 3)
    ├── 06_analyze_box_coupling_and_nulls.py # CENP-B coupling & Stepwise Null Baseline (Fig 4)
    ├── 07_calibrate_length_and_mapping.py   # Physical caliper & MAPQ invariance (Fig 5)
    ├── 08_intra_array_transition.py         # Intra-array epigenetic contrast (Fig 6)
    └── 09_cross_lineage_replication.py      # Cross-lineage biological replication (Fig 7)
```

---

## Reproduction Guide

### Option A: Fast One-Command Reproduction (~15 seconds)
To verify all dynamic ledger assertions and regenerate all 7 publication figures from the authentic data tables:

```bash
# Clone the repository
git clone https://github.com/ad3002/cenpa_nucleosome_architecture.git
cd cenpa_nucleosome_architecture

# Run the master reproduction script
./run_reproduction.sh
```

### Option B: Reviewer Package Standalone Verification
```bash
cd reviewer_package
./reproduce.sh
```

### Option C: Full Raw Pipeline from FASTQ (~1–2 hours)
To download raw FASTQ slices directly from NCBI SRA / EBI European Nucleotide Archive (`SRR13278681`, `SRR13278683`, `SRR13278684`, `SRR15395857`, `SRR9201843`, `SRR9201844`), align with BWA-MEM against CHM13 alpha arrays, and recompute all tables from scratch:

```bash
./run_reproduction.sh --full-raw
```

---

## System Requirements

- **Operating System:** Linux / macOS (tested on Ubuntu 22.04 LTS and macOS Darwin)
- **Software Dependencies:**
  - Python $\ge 3.9$ (`numpy`, `matplotlib`, `scipy`, `openpyxl`)
  - `samtools` $\ge 1.15$
  - `bwa` $\ge 0.7.17$ (for `--full-raw` mode)
  - `pandoc` and `xelatex` (for building PDF manuscripts via `scripts/build_reviewer_package.py`)
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
