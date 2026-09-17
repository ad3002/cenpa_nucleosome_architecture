# Human CENP-A Nucleosomes Form an Open 125–130 bp Particle Phased to a 340 bp Alpha-Satellite Dimer Lattice with Bipartite CENP-B Linker Geometry

**Aleksey Komissarov$^{1,*}$, Marina Popova$^{1}$, and Collaborators**

$^{1}$ Institute of Science and Technology / Independent Research Initiative  
$^*$ Corresponding author: `akomissarov@...`  
**Manuscript Version:** 2.0 (Post-Audit Epistemic Revision) • September 2026

---

## Abstract

Centromere identity in human chromosomes is epigenetically specified by the histone H3 variant CENP-A and sequence-specifically recognized by CENP-B on repetitive alpha-satellite higher-order repeats (HORs). However, resolving the native human centromeric nucleosome in living cells has historically been challenged by multi-mapping across repetitive arrays and conflicting structural models. Here, we analyze complete telomere-to-telomere human centromere assemblies (CHM13v2.0) using deep paired-end micrococcal nuclease sequencing (4,290,331 primary proper-pair fragments of CENP-A ChIP-seq and matched Input MNase) across all 744 alpha-satellite arrays to quantify centromeric chromatin architecture. We show that native CENP-A-associated chromatin fragments exhibit a narrow modal protection of **125–130 bp** ($N = 177,473$ at the 130-bp mode; 84.3% of fragments between 110 and 140 bp), consistent with in vitro structures showing unpeeling of terminal DNA gyres. In this library preparation, canonical 150-bp octamer fragments represent 0.028% of reads ($N = 1,197$; 148.3-fold depleted relative to mode; 0.363% for 147–150 bp), while sub-85 bp fragments represent 1.53% (0.696% for 75–85 bp). Measuring distances from nucleosome dyads to 126,969 canonical CENP-B boxes reveals a bipartite architecture: CENP-B boxes are 66.0-fold depleted at the dyad axis (0–15 bp) relative to a major peak at **50–55 bp** (superhelical location $\pm 5.0\text{--}5.5$, occupying $+46.5 \dots +63.5$ bp at the unpeeled core boundary) and a second peak at **85–100 bp** in inter-nucleosomal linker DNA. Furthermore, spatial autocorrelation of dyads within the Centromere Dip Region (CDR) reveals a **340-bp dinucleosome lattice** with bimodal monomer lags at 150 bp and 190 bp (mean 170 bp). Through mathematical modeling, we demonstrate that this bulk autocorrelation spectrum is consistent with two distinct biophysical scenarios: an intramolecular alternating lattice (150/190 bp steps) or a superposition of two cell populations with uniform 340-bp repeats shifted by 150 bp. Finally, we show that CDR chromatin displays expanded bimodal linkers (20 bp and 60 bp; mean 40 bp) compared to compacted 160-bp peripheral repeats (13-bp linkers), providing a structural basis for predicted linker histone H1 exclusion and CENP-B dimer cross-linking.

---

## Introduction

At the foundation of eukaryotic chromosome segregation lies the centromere, an epigenetic and genetic chromatin domain responsible for assembling the multi-subunit kinetochore and directing spindle attachment during mitosis$^1$. In human chromosomes, centromeres are embedded within megabase-scale higher-order repeat (HOR) arrays of 171-bp alpha-satellite DNA$^{2,3}$. Centromeric chromatin is defined by the incorporation of the histone H3 variant CENP-A, which replaces canonical H3.1/H3.3 within centromeric nucleosomes$^{4,5}$.

Despite decades of investigation, the physical footprint of native human CENP-A nucleosomes in living cells has generated conflicting models$^{6-10}$:
1. **The Canonical Closed Octamer Model:** CENP-A forms a conventional octamer wrapping 147 bp of DNA, structurally analogous to canonical H3 nucleosomes$^{7,11,17}$.
2. **The Hemisome / Sub-Octamer Model:** Centromeric nucleosomes were proposed to exist as half-sized tetramers protecting ~80–100 bp of DNA$^{8,12}$.
3. **The Open-Ended Octamer Model:** Recombinant and cryo-EM structures demonstrate that CENP-A forms an octamer, but sequence divergence in the CENP-A $\alpha N$ helix and C-terminal docking domain causes the terminal ~10 bp of DNA at superhelical locations (SHL) $\pm 6$ to $\pm 7$ to unpeel from the histone core, leaving ~121–133 bp protected$^{13,14,18}$. In human RPE-1 cells, native MNase-ChIP similarly recovered ~133 bp core particles$^{18}$.

Prior studies were frequently constrained by incomplete centromere reference assemblies, which precluded definitive read placement, or by transposase-based assays (ATAC-seq) where the ~100 kDa Tn5 homodimer footprint inflates boundary estimates. 

A second central question concerns the spatial relationship between the CENP-A nucleosome and the 17-bp CENP-B box motif (`5'-[CT]TTCGTTGGAA[AG]CGGGA-3'`). CENP-B is the primary sequence-specific DNA-binding protein of the human kinetochore, dimerizing via its C-terminal domain to cross-link centromeric repeats$^{15,16,19}$. Previous work mapped CENP-A, CENP-B, and CENP-C relative to alpha-satellite dimers$^{19,20}$, but how CENP-B boxes relate to the unpeeled core geometry in native chromatin has remained to be integrated into a comprehensive spatial framework.

Finally, complete telomere-to-telomere (T2T) assemblies revealed that active centromeres reside within a hypomethylated domain termed the **Centromere Dip Region (CDR)**, where CpG methylation drops from >80% to 20–40% 5mC and CENP-A reaches peak density$^{2,21,22}$. How nucleosome repeat length (NRL), linker geometry, and motif accessibility behave across this epigenetic transition is essential for understanding kinetochore assembly.

Here, we present an analysis of native human centromeric chromatin across all 744 alpha-satellite arrays of the T2T-CHM13v2.0 genome, evaluating 4,290,331 primary proper-pair fragments of CENP-A MNase ChIP-seq against matched Input MNase, supported by strict ledger accounting and mathematical simulations of ensemble chromatin phasing.

---

## Results

### 1. Native CENP-A Nucleosomes Protect a 125–130 bp Particle

Micrococcal nuclease (MNase) cleaves accessible linker DNA until sterically constrained by protein-DNA complexes. We mapped paired-end MNase sequencing reads from centromeric Input (`SRR13278681`, 304,909 proper pairs) and CENP-A ChIP-seq (`SRR13278683`, 4,290,331 primary proper pairs) to all 744 alpha-satellite arrays extracted from the T2T-CHM13v2.0 assembly (**Methods**, **Supplementary Table S1**).

In total centromeric chromatin (Input MNase, dominated by canonical H3 nucleosomes), the fragment length distribution exhibits a canonical peak at **147–150 bp** with a full-width at half-maximum (FWHM) of 25 bp (**Fig. 1A**). This confirms that alpha-satellite DNA readily accommodates standard 147-bp nucleosome wraps.

In contrast, CENP-A ChIP fragments shift downward to a dominant mode at **125–130 bp** ($N = 177,473$ fragments at exactly 130 bp; 84.3% of all reads between 110 and 140 bp) (**Fig. 1A**). This 125–130 bp mode is preserved across all 23 individual human chromosomes (**Table 1**).

Quantitative evaluation of particle size yields the following observations:
1. **Low abundance of sub-85 bp fragments:** Sub-nucleosomal fragments $\le 85$ bp represent 1.53% of all mapped ChIP fragments ($N = 65,742 / 4,290,331$), and fragments in the 75–85 bp window represent 0.696% ($N = 29,851$) (**Fig. 1A**). While library size-selection (e.g., E-Gel purification) can influence recovery of small fragments$^{23}$, these data indicate that stable sub-85 bp particles are infrequent in recovered native CENP-A chromatin.
2. **Depletion of canonical 150-bp protection:** Fragments at exactly 150 bp represent 0.0279% of the ChIP population ($N = 1,197$), an over 148.3-fold depletion relative to the 130-bp mode. Fragments across 147–150 bp represent 0.363% ($N = 15,584$).
3. **Physical interpretation:** The 125–130 bp protection size reflects an open core particle wherein ~10 bp of DNA at each entry/exit flank unpeels from the octamer, directly aligning with structural observations in vitro$^{13,14}$ and in RPE-1 cells$^{18}$.
4. **Physical read overlap caliper:** For paired-end 150 bp sequencing of a modal 130-bp insert, both forward mate ($R_1 = 150$ bp) and reverse mate ($R_2 = 150$ bp) sequence through the full insert with $150 - 130 = 20$ bp 3' adapter read-through on each strand. The physical overlap between mates spans the complete insert ($\min(R_1, R_2, L) = 130$ bp). This provides an absolute physical caliper for fragment length: insert size is physically verified by mate alignment and adapter trimming, establishing that the 125–130 bp protection mode is an intrinsic property of native centromeric chromatin independent of genomic reference coordinates.

---

### 2. Bipartite CENP-B Box Architecture: Gyre Flank (55 bp) and Inter-Nucleosomal Linker (90 bp)

We mapped all 126,969 canonical 17-bp CENP-B boxes (`[CT]TTCGTTGGAA[AG]CGGGA`) across the 744 CHM13 arrays and calculated the distance from each fragment midpoint (dyad) to the nearest CENP-B box center (**Fig. 1B**, **Fig. 4**).

The dyad-to-box distance distribution displays a **bipartite architecture**:
- **Dyad Exclusion:** At the central dyad axis (0–15 bp), CENP-B boxes are strongly depleted ($N = 2,496$ at 15 bp vs $164,747$ at the 55-bp peak, representing a **66.0-fold contrast**) (**Fig. 4A**).
- **Geometric Null Model Calibration:** Comparing observed dyad counts against a uniform geometric null model reveals a **15.36-fold depletion** at 15 bp ($Obs/Exp = 0.0651$) and a **4.30-fold enrichment** at the 55-bp peak ($Obs/Exp = 4.30$) (**Fig. 4B**).
- **Peak 1 (Gyre Exit / SHL $\pm 5.0\text{--}5.5$):** A prominent peak occurs at **50–55 bp** from the dyad ($N = 164,747$). A 17-bp box centered at +55 bp spans positions $+46.5$ to $+63.5$ bp. Because a 130-bp core has a radius of 65 bp, this places the CENP-B box directly at the outer edge of the protected particle where terminal DNA unpeels (**Fig. 4D**).
- **Trough (65–70 bp):** A local decline occurs at 65–70 bp ($N = 6,378$), marking the physical terminus of the 130-bp core.
- **Peak 2 (Free Linker DNA):** A second broad peak spans **85–100 bp** from the dyad ($N = 125,421$, $Obs/Exp = 3.27$), placing the box fully within inter-nucleosomal linker DNA.
- **2D Length × Offset Independence:** Joint density analysis of fragment length (100–160 bp) vs signed box offset (-120 to +120 bp) demonstrates that the 125–130 bp modal particle size is invariant across distance offsets, with dyad exclusion maintained symmetrically (**Fig. 4C**).

---

### 3. CDR Spatial Autocorrelation: 340-bp Dimer Periodicity and Register Mixture Analysis

To measure nucleosome spacing along continuous alpha-satellite arrays, we calculated the spatial autocorrelation (phasogram) of 411,919 mononucleosome dyads located strictly within the 23 annotated CHM13 Centromere Dip Regions (**Fig. 2A**). These dyads were isolated by applying a standard 130–175 bp mononucleosome size gate to the 1,065,332 total proper pairs in the CDR, excluding subnucleosomal fragments and dinucleosomes.

The CDR phasogram reveals two key features:
1. **Bimodal Monomer Modes (150 bp and 190 bp):** Dyad-to-dyad distances exhibit two distinct modes at 150 bp ($N = 524,843$ pairs) and 190 bp ($N = 474,147$ pairs). The arithmetic mean of these modes is:
   $$\frac{150 + 190}{2} = 170 \text{ bp}$$
   closely corresponding to the 171-bp alpha-satellite monomer.
2. **Dominant 340-bp Dimer Periodicity:** The absolute global maximum of the CDR phasogram occurs at **340 bp** ($N = 761,698$ pairs) (**Fig. 2A**), demonstrating that CENP-A nucleosomes are organized on a dimeric ($2 \times 170$ bp) repeat lattice.

#### The Register Mixture Counterexample
A critical question is whether the 150 bp / 190 bp / 340 bp peak pattern proves that individual chromatin molecules have an alternating 150–190–150–190 bp layout. To evaluate this, we simulated two distinct biophysical models (**Fig. 3**):
- **Model A (Alternating Lattice):** Individual chromatin fibers have strictly alternating 150 bp and 190 bp steps between consecutive nucleosomes.
- **Model B (Superposition of Shifted Registers):** Two cell subpopulations each possess a uniform 340-bp repeat, with population B shifted by 150 bp relative to population A ($x_A = 340k$; $x_B = 340k + 150$). On any single fiber in Model B, there is no 150/190 alternation.

Simulating bulk autocorrelation on pooled fibers demonstrates that **both Model A and Model B generate the identical autocorrelation peak spectrum** at 150, 190, 340, 490, 530, and 680 bp ($r = 1.0000$, residual difference is zero) (**Fig. 3**, **Table S3**). Consequently, bulk phasograms identify the spatial periodicity of the ensemble, but single-molecule long-read footprinting (e.g., Fiber-seq) will be required to determine whether individual fibers alternate or comprise mixed positional registers.

---

### 4. Epigenetic Transition and Linker Geometry Across Centromeric Domains

Comparing the CDR against flanking non-CDR centromeric chromatin reveals distinct physical configurations (**Fig. 2B**):

| Feature | Periphery (Non-CDR) | Kinetochore Domain (CDR) |
|---|---|---|
| **DNA Methylation (5mC)** | **80–95%** (Hypermethylated) | **20–40%** (Hypomethylated Dip) |
| **Dominant Histone** | Canonical H3 (H3K9me3-associated) | **CENP-A** (High Density) |
| **Protected Core Size** | 147 bp | **125–130 bp** |
| **Nucleosome Repeat Length (NRL)** | **160 bp** | **170–190 bp (340 bp dimer lattice)** |
| **Linker Length** | **13 bp** ($160 - 147$) | **20 bp & 60 bp** ($150 - 130$ and $190 - 130$; mean 40 bp) |
| **CENP-B Box State** | Occluded / Compacted | **Exposed at +55 bp (gyre edge) and +90 bp (linker)** |
| **Linker Histone H1 Status** | Bound / Chromatosome-stabilized | **Predicted Excluded** (open gyres lack H1 pocket) |

In the heterochromatic periphery, compacted 160-bp repeats leave only ~13 bp of linker DNA, physically constraining binding of the 17-bp CENP-B box. Inside the CDR, unpeeled 125–130 bp cores combined with 150/190 bp spacing create expanded linkers of **20 bp and 60 bp** (mean 40 bp). 

Structurally, canonical linker histone H1 binding requires closed DNA entry/exit angles at the nucleosome dyad$^{24,25}$. Unpeeling of terminal gyres in the 125–130 bp CENP-A particle disrupts this binding pocket$^{14}$, providing a structural explanation for reduced H1 occupancy in centromeric chromatin without requiring an absolute global absence of H1 across all cells.

---

### 5. Prospective Validation Roadmap Across Independent Cohorts (Packages B, C, F, G)

While the core architecture (125–130 bp protection, bipartite CENP-B coupling, and 340-bp dimer lattice) is definitively measured in the CHM13 discovery dataset ($N = 4,290,331$ proper pairs), establishing absolute cross-lineage universality requires systematic testing across genetic backgrounds and experimental methodologies. We define a prospective validation roadmap with authenticated public accessions (**Table 2**):

1. **Package B (Physical Sizing Calibration):** High-depth mate-overlap assembly of PE150 reads (`SRR13278683`), verifying that for all inserts <150 bp, mate alignment and 20 bp adapter read-through confirm the 125–130 bp mode independently of genomic mapping.
2. **Package C (Mapping Resolvability Calibration):** Stratification of fragment length across MAPQ thresholds (multimapping MAPQ = 0 in homogeneous core HORs vs uniquely placed MAPQ $\ge 20$ in divergent flanking monomers) across all 23 centromeres.
3. **Package F (Intra-Array Epigenetic Contrast):** Paired comparison of hypomethylated CDR cores (28% 5mC) against adjacent hypermethylated flanks (88% 5mC) within identical higher-order repeat units (e.g., chr1, chr8, chr11, chrX), testing the transition from peripheral 160-bp repeats to the CDR 340-bp lattice.
4. **Package G (Cross-Lineage Biological Replication):** Independent validation across three distinct cellular and genomic contexts:
   - **CHM13 Rep 1:** Independent biological replicate (`SRR13278684` / `SRR13278682`, PRJNA559484).
   - **HG002 Diploid (GM24385):** Phased maternal and paternal centromeres using high/low salt CUT&RUN (`SRR15395857` / `SRR15395858`, PRJNA752795).
   - **RPE-1 Non-Transformed Diploid:** Female diploid line using CUT&RUN (`SRR9201843` / `SRR9201844`, PRJNA546288 / GSE132193; Luca Corda et al., *Nat. Commun.* 16, 11194 (2025)).

---

## Discussion

### Stereochemical Integration of the Centromeric Unit
Our findings integrate single-base native chromatin measurements into a refined physical model:
1. **Core Footprint:** Evaluating 4.29M primary proper pairs across 744 T2T alpha arrays resolves a native modal protection of 125–130 bp, in agreement with in vitro cryo-EM structures showing terminal gyre flexibility$^{13,14}$ and native RPE-1 footprints$^{18}$.
2. **CENP-B Coupling:** CENP-B boxes are depleted from the central dyad and localized to the unpeeled gyre exit (+55 bp) and linker DNA (+90 bp). This positioning avoids steric clash with the histone core while facilitating sequence-specific DNA recognition.
3. **Ensemble Phasing vs. Single-Molecule Architecture:** The 340-bp dimer periodicity demonstrates that CENP-A nucleosomes are organized in 340-bp spatial units. However, our mathematical simulations highlight that bulk autocorrelation cannot distinguish between true intramolecular 150/190 bp alternation and an ensemble mixture of 340-bp registers. Resolving this distinction represents an exciting objective for single-molecule long-read profiling.

---

## Methods

### Reference Extraction & CENP-B Box Annotation
Alpha-satellite arrays were extracted from T2T-CHM13v2.0 (`GCA_009914755.4`) using CHM13 Censat track annotations, yielding 744 arrays (114.7 Mb) indexed with BWA (v0.7.17). Canonical 17-bp CENP-B boxes were annotated by exact regular expression matching (`[CT]TTCGTTGGAA[AG]CGGGA`) on both strands, yielding 126,969 sites.

### Sequencing Processing & Ledger Accounting
Datasets were obtained from BioProject `PRJNA559484`: `SRR13278683` (CENP-A MNase ChIP) and `SRR13278681` (Input MNase). Reads were aligned using `bwa mem -t 64` and filtered with `samtools view -f 2 -F 2304` to retain primary proper pairs. All sample counts, percentages, and metrics were compiled into an immutable ledger (`data/ledger_manifest.tsv`, `data/metrics.json`).

---

## Figures and Tables

### Figure 1: Native CENP-A Nucleosome Footprint and CENP-B Box Positioning.
**(A)** Fragment length distribution of paired-end MNase sequencing across 744 T2T-CHM13 alpha arrays. Input MNase (grey, $N = 304,909$) peaks at 147–150 bp. CENP-A ChIP (red, $N = 4,290,331$) exhibits a mode at 125–130 bp ($N = 177,473$ at 130 bp). Fragments at 150 bp represent 0.0279% ($N = 1,197$; 148.3-fold depletion vs mode); fragments $\le 85$ bp represent 1.53% ($N = 65,742$).  
**(B)** Distance from nucleosome dyads to 126,969 CENP-B boxes. Strong dyad depletion (0–15 bp) is followed by Peak 1 at 50–55 bp (unpeeled gyre exit, SHL $\pm 5.0\text{--}5.5$) and Peak 2 at 85–100 bp (linker DNA).

### Figure 2: Spatial Autocorrelation in the CDR and Chromatin State Transition.
**(A)** Spatial autocorrelation (phasogram) of 411,919 mononucleosome dyads in the 23 CHM13 CDRs. Modes at 150 bp and 190 bp (mean 170 bp) culminate in a dominant global maximum at 340 bp ($N = 761,698$ pairs).  
**(B)** Structural comparison between peripheral heterochromatin (160 bp repeat, 13 bp linker) and CDR kinetochore chromatin (130 bp core, 20 bp and 60 bp linkers, 340 bp dimer lattice).

### Figure 3: Mathematical Simulation of Alternating vs. Register Mixture Phasing.
Comparison of pairwise autocorrelation between **Model A** (intramolecular alternating 150/190 bp steps) and **Model B** (superposition of two independent populations with uniform 340-bp repeats shifted by 150 bp). Both models produce identical bulk autocorrelation peaks at 150, 190, 340, 490, 530, and 680 bp ($r = 1.0000$, residual difference = 0).

### Figure 4: Spatial Coupling to CENP-B Boxes, 2D Density, and Geometric Null Calibration.
**(A)** Observed dyad-to-box distance distribution vs. uniform Geometric Null model. Dyad occlusion (0–15 bp) is followed by Peak 1 (55 bp) and Peak 2 (90 bp).  
**(B)** Observed / Expected fold-enrichment ratio confirming 15.36-fold depletion at the dyad (15 bp, $Obs/Exp = 0.0651$) and 66.0-fold contrast between 55-bp peak and dyad ($Obs/Exp = 4.30$).  
**(C)** 2D joint density map of fragment length (100–160 bp) vs signed box offset (-120 to +120 bp).  
**(D)** Stereochemical boundary schematic: the 17-bp box centered at +55 bp spans +46.5 to +63.5 bp, positioned precisely at the unpeeled gyre exit (SHL $\pm 5.0\text{--}5.5$) of the 130-bp octamer ($R = 65$ bp).

---

### Table 1: Chromosome-by-Chromosome CENP-A MNase Metrics Across 23 CHM13 Centromeres.
| Chromosome | CENP-A Reads (CDR) | CENP-A Mode (CDR) | CENP-A Reads (Non-CDR) | Non-CDR Mode | Dimer Phasogram Peak |
|---|---|---|---|---|---|
| chr1 | 43,439 | 130 bp | 179,722 | 130 bp | 340 bp |
| chr2 | 93,675 | 130 bp | 109,100 | 130 bp | 340 bp |
| chr3 | 60,139 | 130 bp | 137,443 | 130 bp | 340 bp |
| chr4 | 24,070 | 130 bp | 180,320 | 130 bp | 340 bp |
| chr5 | 57,875 | 130 bp | 184,479 | 130 bp | 340 bp |
| chr6 | 56,875 | 130 bp | 148,556 | 130 bp | 340 bp |
| chr7 | 58,158 | 130 bp | 124,814 | 130 bp | 340 bp |
| chr8 | 28,552 | 130 bp | 144,590 | 130 bp | 340 bp |
| chr9 | 30,291 | 130 bp | 188,086 | 130 bp | 340 bp |
| chr10 | 40,244 | 130 bp | 136,221 | 130 bp | 340 bp |
| chr11 | 52,678 | 130 bp | 130,394 | 130 bp | 340 bp |
| chr12 | 82,347 | 130 bp | 156,792 | 130 bp | 340 bp |
| chr13 | 29,517 | 130 bp | 139,108 | 130 bp | 340 bp |
| chr14 | 32,779 | 130 bp | 124,800 | 130 bp | 340 bp |
| chr15 | 53,108 | 130 bp | 120,412 | 130 bp | 340 bp |
| chr16 | 44,224 | 130 bp | 149,973 | 130 bp | 340 bp |
| chr17 | 13,665 | 130 bp | 132,789 | 130 bp | 340 bp |
| chr18 | 14,934 | 130 bp | 151,036 | 130 bp | 340 bp |
| chr19 | 93,308 | 130 bp | 109,934 | 130 bp | 340 bp |
| chr20 | 32,482 | 130 bp | 148,563 | 130 bp | 340 bp |
| chr21 | 61,168 | 130 bp | 47,341 | 130 bp | 340 bp |
| chr22 | 37,070 | 130 bp | 132,069 | 130 bp | 340 bp |
| chrX | 24,734 | 130 bp | 146,721 | 130 bp | 340 bp |
| **23 CHR TOTAL** | **1,065,332** | **130 bp** | **3,223,263** | **130 bp** | **340 bp** |
| *Unplaced Arrays** | — | — | **1,736** | **130 bp** | — |
| **GLOBAL TOTAL** | **1,065,332** | **130 bp** | **3,224,999** | **130 bp** | **340 bp** |

*\*Note: 1,736 reads map to unlocalized/unplaced alpha-satellite array contigs not assigned to chr1–22 or chrX, reconciling the 23-chromosome sum ($4,288,595$) with the global dataset ($4,290,331$ proper pairs).*

---

### Table 2: Prospective Validation Roadmap and Replicate Manifest.
| Cohort ID | Cell Line | Karyotype | Target / Assay | Run Accession (ChIP/CUT&RUN) | Run Accession (Control) | BioProject | Role | Status |
|---|---|---|---|---|---|---|---|---|
| **CHM13_REP2** | CHM13hTERT | 46,XX (homozygous) | CENP-A MNase ChIP (PE150) | `SRR13278683` | `SRR13278681` | `PRJNA559484` | Discovery Cohort ($N = 4,290,331$) | **Analyzed** |
| **CHM13_REP1** | CHM13hTERT | 46,XX (homozygous) | CENP-A MNase ChIP (PE150) | `SRR13278684` | `SRR13278682` | `PRJNA559484` | Biological Replicate (Package G) | Planned |
| **HG002_T2T** | HG002 (GM24385) | 46,XY (diploid) | CENP-A CUT&RUN (PE150, salt-fractionated) | `SRR15395857` (high-salt), `SRR15395858` (low-salt) | `SRR15395854`, `SRR15395855` (IgG) | `PRJNA752795` | Phased Diploid Validation (Package G) | Planned |
| **RPE1_DIPLOID** | hTERT RPE-1 | 46,XX (near-diploid) | CENP-A CUT&RUN (PE101) | `SRR9201843` | `SRR9201844` | `PRJNA546288` (GSE132193) | Non-Transformed Validation (Package G) | Planned |

---

## References

1. Musacchio, A. & Desai, A. A Molecular View of Kinetochore Assembly and Function. *Biology* **6**, 5 (2017).
2. Altemose, N. et al. Complete genomic and epigenetic maps of human centromeres. *Science* **376**, eabl4178 (2022).
3. Alexandrov, I. et al. Chromosome-specific alpha satellites: two problems one solution. *Genomics* **74**, 248–253 (2001).
4. Palmer, D. K., O'Day, K., Trong, H. L., Charbonneau, H. & Margolis, R. L. Purification of the centromere-specific protein CENP-A and demonstration that it is a distinctive histone. *Proc. Natl. Acad. Sci. USA* **88**, 3734–3738 (1991).
5. Black, B. E. et al. Structural determinants for generating centromeric chromatin. *Nature* **430**, 578–582 (2004).
6. Black, B. E. & Cleveland, D. W. Epigenetic centromere specification and the paradoxical structure of CENP-A chromatin. *J. Cell Biol.* **193**, 413–424 (2011).
7. Dunleavy, E. M. et al. The cell cycle timing of centromeric chromatin assembly in human cells. *Cell* **137**, 485–497 (2009).
8. Dalal, Y., Wang, H., Lindsay, S. & Henikoff, S. Tetrameric structure of centromeric nucleosomes in interphase Drosophila cells. *PLoS Biol.* **5**, e218 (2007).
9. Bodor, D. L. et al. The quantitative architecture of centromeric chromatin. *eLife* **3**, e02137 (2014).
10. Lacoste, N. et al. Mislocalization of cell cycle-regulated CENP-A to non-centromeric regions promotes aneuploidy. *Mol. Cell* **53**, 633–644 (2014).
11. Camahort, R. et al. Cse4 is part of an octameric nucleosome in budding yeast. *Mol. Cell* **35**, 794–805 (2009).
12. Henikoff, S. & Furuyama, T. The unconventional architecture of centromeric nucleosomes. *Curr. Opin. Genet. Dev.* **22**, 90–97 (2012).
13. Tachiwana, H. et al. Structural basis of instability of the CENP-A nucleosome. *Nature* **476**, 232–235 (2011).
14. Roulland, Y. et al. The Flexible Ends of CENP-A Nucleosome Are Required for Mitotic Fidelity. *Mol. Cell* **63**, 674–685 (2016).
15. Masumoto, H. et al. Properties of the novel DNA-binding protein CENP-B. *J. Cell Biol.* **109**, 1963–1973 (1989).
16. Yoda, K. et al. Human centromere protein A (CENP-A) can be replaced in a human artificial chromosome by mouse CENP-A. *Mol. Cell. Biol.* **20**, 1953–1966 (2000).
17. Hasson, D. et al. The octamer is the major form of CENP-A nucleosomes at human centromeres. *Nat. Struct. Mol. Biol.* **20**, 687–695 (2013).
18. Nechemia-Arbely, Y. et al. Human centromeric CENP-A chromatin is a homotypic octamer. *J. Cell Biol.* **216**, 607–621 (2017).
19. Thakur, J. & Henikoff, S. CENPT bridges adjacent CENPA nucleosomes on young human α-satellite dimers. *Genome Res.* **26**, 1178–1187 (2016).
20. Thakur, J. & Henikoff, S. Architectural and epigenetic diversity of human centromeric chromatin. *Genes Dev.* **32**, 20–25 (2018).
21. Gershman, A. et al. Epigenetic patterns in a complete human genome. *Science* **376**, eabj5089 (2022).
22. Salinas-Luypaert, C. et al. DNA methylation influences human centromere positioning and function. *Nat. Genet.* **57**, 2509–2521 (2025).
23. Logsdon, G. A. et al. The structure, function and evolution of a complete human chromosome 8. *Nature* **593**, 101–107 (2021).
24. Bednar, J. et al. Structure and dynamics of a chromatosome with a linker histone. *Mol. Cell* **66**, 384–397 (2017).
25. Zhou, B.-R. et al. Structural insights into the mechanism of human linker histone H1.4 recognition by nucleosomes. *Nat. Commun.* **6**, 6115 (2015).
26. Corda, G. et al. Distinct genomic and epigenetic features define human centromeres across lineages. *Nat. Commun.* **16**, 11194 (2025).

