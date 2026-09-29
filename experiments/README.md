# CENP-A Nucleosome Architecture — Dedicated Experimental Suite (EXP)

This directory contains standalone, reproducible computational experiments designed to test mechanistic predictions of the CENP-A chromatin architecture model against independent datasets and physical interventions.

Symlinked as `EXP/` at the repository root.

---

## Experimental Registry

| Experiment ID | Directory | Description | Datasets | Status | Key Finding |
| :--- | :--- | :--- | :--- | :---: | :--- |
| **EXP-03** | [`exp03_dyad_vs_ends_geometry/`](./exp03_dyad_vs_ends_geometry/) | Resolving the geometric anchor of CENP-B coupling: dyad center vs fragment ends | CHM13 Rep 2 (`SRR13278683`), Rep 1 (`SRR13278684`) | **COMPLETE** | **H₁ confirmed:** Central dyad is fixed ($\beta \approx +0.10$); terminal arms unpeel symmetrically ($-0.39$ and $+0.61$); H₂ ($\beta = \pm 0.5$) and H₃ (box switching) falsified ($p < 10^{-15}$). |
| **EXP-04** | [`exp04_cenpb_sequence_controls/`](./exp04_cenpb_sequence_controls/) | Specificity of CENP-B coupling under sequence background & cleavage bias controls ($M_0$ vs $M_1$) | CHM13 CENP-A/Input (`SRR13278681`), RPE-1 CENP-B (`SRR9201844`) | **COMPLETE** | $M_1$ decisively preferred over $M_0$ ($p < 10^{-15}$, $\Delta\text{BIC} > 320,000$); Input MNase lacks +55/+95 bp peaks; monotonic dose-response attenuation across $d_H=0 \to 1 \to 2 \to B^-$; direct RPE-1 CENP-B footprint confirmed. |
| **EXP-02** | [`exp02_single_molecule_fiberseq/`](./exp02_single_molecule_fiberseq/) | Single-molecule Fiber-seq testing of nucleosome spacing alternation ($P(g_i, g_{i+1})$, $\text{corr}(g_i, g_{i+1})$) | CHM13 `GSM7074431` (`GSE226394`, 1.34M pairs) | **COMPLETE** | Resolves Phasogram Degeneracy: **$H_2$ Confirmed (Register Mixture Model)** ($r = +0.0781, p < 10^{-136}$; variance ratio $= 1.074 > 1$); $H_1$ (alternating lattice, $r < 0$) refuted; open $115-135$ bp core enriched in CDR ($21.7\%$ vs $16.3\%$ canonical). |
| **EXP-01** | [`exp01_demethylation_dynamics/`](./exp01_demethylation_dynamics/) | Chromatin geometry response to targeted demethylation | Salinas-Luypaert 2025 (`PRJNA1270043` / Zenodo 15875037) | **COMPLETE** | **$H_1$ Confirmed:** Local core geometry is strictly invariant ($128$ bp mode, $\Delta L = 0.18$ bp $< 2.0$ bp); $H_2$ falsified; **$H_3$ Confirmed:** linker accessibility increases significantly ($25.1\% \to 38.2\%$, $p < 10^{-15}$); $+45.2$ kb outward domain expansion past CDR boundary. |
| **EXP-05** | [`exp05_solubilization_salt_titration/`](./exp05_solubilization_salt_titration/) | Disentangling stable short core from extraction fractionation bias | Thakur & Henikoff (`GSE104805`, 0–500 mM) & HG002 | **COMPLETE** | **$H_1$ Confirmed:** 127–128 bp open core persists across independent soluble chemistries (0 mM ChIP & high-salt CUT&RUN); $H_2$ falsified (not a 0 mM artifact); **$H_3$ Confirmed:** central dyad phasing is strictly invariant ($r = 0.672, p < 10^{-5}$) across extraction conditions while less soluble fractions capture extended terminal protection. |
| **EXP-06** | [`exp06_structural_accessibility_model/`](./exp06_structural_accessibility_model/) | 3D steric collision modeling of CENP-B accessibility | PDB `1HLV`, `6SE0`, `1KX5`, Nagpal 2023 | **COMPLETE** | **$H_3$ Confirmed:** Conformational unpeeling ensemble $A(d, u, \theta)$ predicts bimodal genomic positioning ($r = 0.376, \Delta\text{BIC} > 10$) without parameter tuning; spontaneous 10–15 bp unpeeling relieves clashes at $+55$ bp exit site; $+95$ bp linker is unconstrained. |

---

## Standard Structure for Each Experiment Module

Every experiment module `expXX_*/` contains:
1. `run_experiment.sh` — 1-command executable runner.
2. `analyze_*.py` — Standalone, self-contained Python analysis script.
3. `data/` — Output TSV and JSON summaries with bootstrap confidence intervals.
4. `figures/` — Publication-ready figures (PNG, SVG, PDF).
5. `REPORT.md` — Formal scientific report detailing objectives, competing hypotheses, empirical results, and hypothesis decisions.
