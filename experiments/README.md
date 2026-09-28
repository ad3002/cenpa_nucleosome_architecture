# CENP-A Nucleosome Architecture — Dedicated Experimental Suite (EXP)

This directory contains standalone, reproducible computational experiments designed to test mechanistic predictions of the CENP-A chromatin architecture model against independent datasets and physical interventions.

Symlinked as `EXP/` at the repository root.

---

## Experimental Registry

| Experiment ID | Directory | Description | Datasets | Status | Key Finding |
| :--- | :--- | :--- | :--- | :---: | :--- |
| **EXP-03** | [`exp03_dyad_vs_ends_geometry/`](./exp03_dyad_vs_ends_geometry/) | Resolving the geometric anchor of CENP-B coupling: dyad center vs fragment ends | CHM13 Rep 2 (`SRR13278683`), Rep 1 (`SRR13278684`) | **COMPLETE** | **H₁ confirmed:** Central dyad is fixed ($\beta \approx +0.10$); terminal arms unpeel symmetrically ($-0.39$ and $+0.61$); H₂ ($\beta = \pm 0.5$) and H₃ (box switching) falsified ($p < 10^{-15}$). |
| **EXP-04** | [`exp04_cenpb_sequence_controls/`](./exp04_cenpb_sequence_controls/) | Specificity of CENP-B coupling under sequence background & cleavage bias controls ($M_0$ vs $M_1$) | CHM13 CENP-A/Input (`SRR13278681`), RPE-1 CENP-B (`SRR9201844`) | **COMPLETE** | $M_1$ decisively preferred over $M_0$ ($p < 10^{-15}$, $\Delta\text{BIC} > 320,000$); Input MNase lacks +55/+95 bp peaks; monotonic dose-response attenuation across $d_H=0 \to 1 \to 2 \to B^-$; direct RPE-1 CENP-B footprint confirmed. |
| **EXP-02** | [`exp02_single_molecule_fiberseq/`](./exp02_single_molecule_fiberseq/) | Single-molecule Fiber-seq testing of nucleosome spacing alternation ($P(g_i, g_{i+1})$, $\text{corr}(g_i, g_{i+1})$) | CHM13 `GSM7074431` (`GSE226394`) | *Queued (Wave 2)* | Direct single-molecule test of alternating lattice (H₁) vs register mixture (H₂) vs irregular packing (H₃). |
| **EXP-01** | [`exp01_demethylation_dynamics/`](./exp01_demethylation_dynamics/) | Chromatin geometry response to targeted demethylation | Salinas-Luypaert 2025 (`PRJNA1270043`) | *Queued (Wave 2)* | Tests whether domain expansion preserves local geometry (H₁), shifts parameters (H₂), or alters state weights (H₃). |
| **EXP-05** | [`exp05_solubilization_salt_titration/`](./exp05_solubilization_salt_titration/) | Disentangling stable short core from extraction fractionation bias | Thakur & Henikoff (`GSE104805`, 0–500 mM) | *Queued (Wave 3)* | Tests whether the 133-bp core is an intrinsic feature across salt fractions or a soluble-specific state. |
| **EXP-06** | [`exp06_structural_accessibility_model/`](./exp06_structural_accessibility_model/) | 3D steric collision modeling of CENP-B accessibility | PDB `1HLV`, `6SE0`, `1KX5`, Nagpal 2023 | *Queued (Wave 3)* | Establishes steric collision map $A(d, u, \theta)$ and predicts genomic accessibility rules. |

---

## Standard Structure for Each Experiment Module

Every experiment module `expXX_*/` contains:
1. `run_experiment.sh` — 1-command executable runner.
2. `analyze_*.py` — Standalone, self-contained Python analysis script.
3. `data/` — Output TSV and JSON summaries with bootstrap confidence intervals.
4. `figures/` — Publication-ready figures (PNG, SVG, PDF).
5. `REPORT.md` — Formal scientific report detailing objectives, competing hypotheses, empirical results, and hypothesis decisions.
