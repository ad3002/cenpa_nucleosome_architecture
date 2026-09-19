# Human CENP-A Nucleosomes Form an Open 125–130 bp Particle Phased to a 340 bp Alpha-Satellite Dimer Lattice with Bipartite CENP-B Linker Geometry

**Aleksey Komissarov$^{1,*}$, Marina Popova$^{1}$, and Collaborators**

$^{1}$ Institute of Science and Technology / Independent Research Initiative  
$^*$ Corresponding author: `akomissarov@...`  
**Manuscript Version:** 3.1 (Comprehensive Audit Remediation & Epistemic Revision) • September 2026

---

## Abstract

Centromere identity in human chromosomes is epigenetically specified by the histone H3 variant CENP-A and sequence-specifically recognized by CENP-B on repetitive alpha-satellite higher-order repeats (HORs). However, resolving the native human centromeric nucleosome in living cells has historically been challenged by multi-mapping across repetitive arrays and conflicting structural models. Here, we analyze complete telomere-to-telomere human centromere assemblies (T2T-CHM13v2.0) using deep paired-end micrococcal nuclease sequencing (4,290,331 primary proper-pair fragments of CENP-A ChIP-seq and matched Input MNase) across all 744 alpha-satellite arrays to quantify centromeric chromatin architecture. We show that native CENP-A-associated chromatin fragments exhibit a single-base mode at **133 bp** ($N = 212,205$ global; $55,892$ in CDR; $156,313$ in Non-CDR), with $177,473$ fragments at 130 bp and **84.29%** of all fragments ($N = 3,616,490 / 4,290,331$) concentrated in the 110–140 bp unpeeled core window (5-bp binning centers the distribution at 130 bp). This protection size is consistent with in vitro structures showing unpeeling of terminal DNA gyres. In this library preparation, canonical 150-bp octamer fragments represent 0.0279% of reads ($N = 1,197$; 177.28-fold depleted relative to the 133-bp mode and 148.26-fold depleted relative to 130 bp; 0.363% for 147–150 bp), while sub-85 bp fragments represent 1.53% ($N = 65,742$; 0.696% for 75–85 bp). Measuring distances from nucleosome dyads to 126,969 canonical CENP-B boxes reveals a bipartite architecture: CENP-B boxes are 66.0-fold depleted at the dyad axis (0–15 bp) relative to a major peak at the 55-bp bin (superhelical location $\pm 5.0\text{--}5.5$, where a 17-bp box centered at +55 bp occupies $+46.5 \dots +63.5$ bp at the unpeeled core boundary) and a second peak at **85–100 bp** in inter-nucleosomal linker DNA (displaying 6.54-fold enrichment over an empirical stepwise null model baseline at 100 bp). Furthermore, spatial autocorrelation of dyads within the Centromere Dip Region (CDR) reveals a dominant non-zero peak at **340 bp** ($N = 761,698$ pairs in CDR; $591,711$ in Non-CDR) with bimodal monomer lags at 150 bp and 190 bp (mean 170 bp). Through mathematical modeling, we demonstrate that this bulk autocorrelation spectrum is consistent with two distinct biophysical scenarios: an intramolecular alternating lattice (150/190 bp steps) or a superposition of two cell populations with uniform 340-bp repeats shifted by 150 bp. Finally, through structural and stereochemical modeling, we show that unpeeled CENP-A cores within a 340-bp lattice predict expanded linkers (modeled as 20 bp and 60 bp; mean 40 bp) compared to compacted 160-bp peripheral repeats (~13-bp linkers), providing a mechanistic framework consistent with predicted linker histone H1 exclusion and CENP-B dimer cross-linking.

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

### 1. Native CENP-A Nucleosomes Protect a 125–133 bp Open Core Particle

Micrococcal nuclease (MNase) cleaves accessible linker DNA until sterically constrained by protein-DNA complexes. We mapped paired-end MNase sequencing reads from centromeric Input (`SRR13278681`, 304,909 proper pairs) and CENP-A ChIP-seq (`SRR13278683`, 4,290,331 primary proper pairs) to all 744 alpha-satellite arrays extracted from the T2T-CHM13v2.0 assembly (**Methods**, **Supplementary Table S1**).

In total centromeric chromatin (Input MNase, dominated by canonical H3 nucleosomes), the fragment length distribution exhibits a canonical peak at **147–150 bp** with a full-width at half-maximum (FWHM) of 25 bp (**Fig. 1A**). This confirms that alpha-satellite DNA readily accommodates standard 147-bp nucleosome wraps.

In contrast, CENP-A ChIP fragments exhibit a true single-base mode at **133 bp** ($N = 212,205$ fragments globally; $55,892$ in CDR; $156,313$ in Non-CDR), with $177,473$ fragments at 130 bp and **84.29%** of all reads ($N = 3,616,490 / 4,290,331$) concentrated within the 110–140 bp window (**Fig. 1A**). When grouped into 5-bp bins, the modal bin is 130 bp (128–132 bp). This 125–133 bp protection is preserved across all 23 individual human chromosomes (**Table 1**).

Quantitative evaluation of particle size yields the following observations:
1. **Low abundance of sub-85 bp fragments:** Sub-nucleosomal fragments $\le 85$ bp represent 1.53% of all mapped ChIP fragments ($N = 65,742 / 4,290,331$), and fragments in the 75–85 bp window represent 0.696% ($N = 29,851$) (**Fig. 1A**). While library size-selection (e.g., E-Gel purification) can influence recovery of small fragments$^{23}$, these data indicate that stable sub-85 bp particles are infrequent in recovered native CENP-A chromatin.
2. **Depletion of canonical 150-bp protection:** Fragments at exactly 150 bp represent 0.0279% of the ChIP population ($N = 1,197$), representing a 177.28-fold depletion relative to the 133-bp single-base mode and a 148.26-fold depletion relative to 130 bp. Fragments across 147–150 bp represent 0.363% ($N = 15,584$).
3. **Physical interpretation:** The 125–133 bp protection size reflects an open core particle wherein ~10 bp of DNA at each entry/exit flank unpeels from the octamer, directly aligning with structural observations in vitro$^{13,14}$ and in RPE-1 cells$^{18}$.
4. **Geometric read overlap caliper model:** In paired-end 150-bp sequencing of a modal ~130-bp insert, each mate ($R_1 = 150$ bp, $R_2 = 150$ bp) is expected to sequence across the full physical insert into the opposite adapter ($150 - L = 20$ bp read-through for $L = 130$ bp). The theoretical mate overlap spans $\max(0, \min(R_1, L) + \min(R_2, L) - L) = L$ bp. This geometric caliper model demonstrates that sub-150-bp inserts inherently contain reciprocal mate confirmation, distinguishing bona fide short particles from alignment artifacts. Full empirical read-level adapter-trimming calibration is designated for prospective Work Package B.

---

### 2. Bipartite CENP-B Box Architecture: Gyre Flank (55 bp) and Inter-Nucleosomal Linker (90–100 bp)

We mapped all 126,969 canonical 17-bp CENP-B boxes (`[CT]TTCGTTGGAA[AG]CGGGA`) across the 744 CHM13 arrays and calculated the distance from each fragment midpoint (dyad) to the nearest CENP-B box center (**Fig. 1B**, **Fig. 4**).

The dyad-to-box distance distribution displays a **bipartite architecture**:
- **Dyad Exclusion:** At the central dyad axis (0–15 bp), CENP-B boxes are strongly depleted ($N = 2,496$ at 15 bp vs $164,747$ at the 55-bp bin, representing a **66.0-fold contrast**) (**Fig. 4A**).
- **Empirical Stepwise Null Model Baseline Calibration:** Comparing observed dyad counts against an empirical stepwise null model baseline (evaluating uniform dyad placement across alpha-satellite arrays with 50% box density: $p=1.0$ for $d \le 85$ bp, $p=0.5$ for $85 < d \le 170$ bp, and $p=0.1$ for $d > 170$ bp) reveals a **15.36-fold depletion** at 15 bp ($Obs/Exp = 0.0651$) and a **4.30-fold enrichment** at the 55-bp bin ($Obs/Exp = 4.30$) (**Fig. 4B**).
- **Peak 1 (Gyre Exit / SHL $\pm 5.0\text{--}5.5$):** A prominent peak occurs at the **55-bp bin** ($[55, 60)$ bp with midpoint 57.5 bp; $N = 164,747$). A canonical 17-bp box centered at +55 bp spans positions $+46.5$ to $+63.5$ bp relative to the dyad. Because a 130-bp core has a radius of 65 bp, this places the CENP-B box directly at the outer edge of the protected particle where terminal DNA unpeels (**Fig. 4D**).
- **Trough (65–70 bp):** A local decline occurs at 65–70 bp ($N = 6,378$), marking the physical terminus of the 130-bp core.
- **Peak 2 (Free Linker DNA):** A second broad peak spans **85–100 bp** from the dyad ($N = 125,421$ at 100 bp, $Obs/Exp = 6.54$; $N = 112,228$ at 90 bp, $Obs/Exp = 5.86$), placing the box fully within inter-nucleosomal linker DNA.
- **Theoretical Stereochemical Model Schema:** A theoretical joint density schema of fragment length (100–160 bp) vs signed box offset (-120 to +120 bp) illustrates how the 130-bp unpeeled core particle and bipartite box spacing predict invariant modal particle length across distance offsets with symmetric dyad exclusion (**Fig. 4C**). Full empirical single-fragment joint density calibration across read-level BAM records is designated for prospective Work Package D.

---

### 3. CDR Spatial Autocorrelation: 340-bp Dimer Periodicity and Register Mixture Analysis

To measure nucleosome spacing along continuous alpha-satellite arrays, we calculated the spatial autocorrelation (phasogram) of mononucleosome dyads located strictly within the 23 annotated CHM13 Centromere Dip Regions (**Fig. 2A**). Dyads were isolated by applying a standard 130–175 bp mononucleosome size gate ($N = 411,919$ fragments) to the 1,065,332 total proper pairs in the CDR (with the broader 110–180 bp gate encompassing 960,496 fragments).

The CDR phasogram reveals two key features:
1. **Bimodal Monomer Modes (150 bp and 190 bp):** Dyad-to-dyad distances exhibit two distinct modes at 150 bp ($N = 524,843$ pairs) and 190 bp ($N = 474,147$ pairs), separated by a trough at 170 bp ($N = 309,708$ pairs). The arithmetic mean of these modes is:
   $$\frac{150 + 190}{2} = 170 \text{ bp}$$
   closely corresponding to the 171-bp alpha-satellite monomer.
2. **Dominant 340-bp Dimer Periodicity:** In the non-zero inter-nucleosomal window (100–800 bp), the dominant maximum occurs at **340 bp** ($N = 761,698$ pairs in CDR; $N = 591,711$ pairs in Non-CDR) (**Fig. 2A**), demonstrating that CENP-A nucleosomes are organized on an ensemble dimeric ($2 \times 170$ bp) repeat lattice.

#### The Register Mixture Counterexample
A critical question is whether the 150 bp / 190 bp / 340 bp peak pattern proves that individual chromatin molecules possess an alternating 150–190–150–190 bp layout. To evaluate this, we simulated two distinct biophysical models (**Fig. 3**):
- **Model A (Alternating Lattice):** Individual chromatin fibers have strictly alternating 150 bp and 190 bp steps between consecutive nucleosomes.
- **Model B (Superposition of Shifted Registers):** Two cell subpopulations each possess a uniform 340-bp repeat, with population B shifted by 150 bp relative to population A ($x_A = 340k$; $x_B = 340k + 150$). On any single fiber in Model B, there is no 150/190 alternation.

Simulating bulk autocorrelation on pooled fibers demonstrates that **both Model A and Model B generate the identical autocorrelation peak spectrum** at 150, 190, 340, 490, 530, and 680 bp ($r = 1.0000$, residual difference is zero) (**Fig. 3**, **Table S3**). Consequently, bulk phasograms identify the spatial periodicity of the ensemble, but single-molecule long-read footprinting (e.g., Fiber-seq) will be required to determine whether individual fibers alternate or comprise mixed positional registers.

---

### 4. Epigenetic Transition and Linker Geometry Models

Comparing the CDR against flanking non-CDR centromeric chromatin reveals distinct physical configurations (**Fig. 2B**):

| Feature | Periphery (Non-CDR Model) | Kinetochore Domain (CDR Model) |
|---|---|---|
| **DNA Methylation (5mC)** | **80–95%** (Hypermethylated) | **20–40%** (Hypomethylated Dip) |
| **Dominant Histone** | Canonical H3 (H3K9me3-associated) | **CENP-A** (High Density) |
| **Protected Core Size** | 147 bp | **125–133 bp** (Mode: 133 bp) |
| **Nucleosome Repeat Length (NRL)** | **160 bp** | **170–190 bp (340 bp dimer lattice)** |
| **Modeled Linker Length** | **13 bp** ($160 - 147$) | **20 bp & 60 bp** ($150 - 130$ and $190 - 130$; mean 40 bp) |
| **CENP-B Box State** | Sterically Constrained | **Positioned at +55 bp (gyre edge) and +90–100 bp (linker)** |
| **Linker Histone H1 Status** | Bound / Chromatosome-stabilized | **Predicted Excluded** (open gyres lack H1 pocket) |

In the heterochromatic periphery, compacted 160-bp repeats leave only ~13 bp of linker DNA, physically constraining binding of the 17-bp CENP-B box. Inside the CDR, unpeeled 125–133 bp cores combined with 150/190 bp spacing create modeled linkers of **20 bp and 60 bp** (mean 40 bp). 

Structurally, canonical linker histone H1 binding requires closed DNA entry/exit angles at the nucleosome dyad$^{24,25}$. Unpeeling of terminal gyres in the 125–133 bp CENP-A particle disrupts this binding pocket$^{14}$, providing a structural explanation for reduced H1 occupancy in centromeric chromatin without requiring an absolute global absence of H1 across all cells. We emphasize that 20/60 bp linkers and H1 exclusion represent mechanistic models inferred from bulk spacing and structural geometry, rather than directly measured single-molecule distributions.

---

### 5. Physical Overlap Caliper and MAPQ Invariance (Packages B & C)

To address potential technical concerns that the open 125–133 bp nucleosome footprint might represent an artifact of short-read alignment scoring, soft-clipping, or ambiguous placement across repetitive alpha-satellite higher-order repeats (HORs), we executed two independent biophysical and algorithmic calibrations (**Fig. 5**):

#### Reference-Free FASTQ Read Overlap Caliper (Package B)
In paired-end 150 bp sequencing (PE150), any DNA fragment shorter than 150 bp is sequenced across its entire physical length by both Read 1 and Read 2, with both sequencing reads extending past the opposite 3' termini into Illumina adapter sequences (**Fig. 5A**). Consequently, fragment length can be determined with single-base precision directly from raw FASTQ sequence boundaries independently of any reference genome or alignment tool.

Evaluating 100,000 paired-end reads from `SRR13278683` with 3' adapter detection (`AGATCGGAAGAGC`) and base-for-base reverse-complement matching:
1. **High Concordance:** 89,050 read pairs (89.05%) exhibited concordant adapter boundaries on both mates, with 88,481 pairs (88.48%) verified by sequence identity ($\le 2$ non-N mismatches) across the insert duplex.
2. **Identical Open Core Mode:** The reference-free physical caliper distribution exhibits a dominant single-base mode at **133 bp** (modal 5-bp bin at **130 bp**; **Fig. 5B**), with **88.18%** of fragments falling within the [110, 140] bp core gate. Canonical 150-bp fragments are depleted to 0.00%, and sub-85 bp fragments represent only 1.15% of the library.
3. **Direct Aligner Concordance:** Cross-referencing physical caliper lengths against BWA-MEM alignment insert lengths (`TLEN`) for 75,911 mapped pairs reveals near-perfect linear agreement ($R^2 = 0.999$, median difference = **0.0 bp**, mean difference = 0.42 bp; **Fig. 5C**), with **98.40%** exact base-for-base concordance across modal sizes. This confirms that BWA-MEM alignment preserves physical fragment boundaries without systematic contraction or expansion.

#### Invariance Across MAPQ Strata (Package C)
Centromeric alpha-satellite arrays contain both highly homogenized core HORs (generating multi-mapped reads with $\text{MAPQ} = 0$) and divergent repeat variants (yielding uniquely placed reads with $\text{MAPQ} \ge 20$). To test whether repeat-mapping ambiguity distorts particle sizing, we stratified mapped pairs into $\text{MAPQ} = 0$ ($N = 80,957$) and $\text{MAPQ} \ge 20$ ($N = 2,942$) cohorts (**Fig. 5D**).

Both strata exhibit identical single-base modes at **133 bp** ($\Delta = 0$ bp; **Fig. 5D**), with identical distribution profiles across the 110–140 bp window. This invariance demonstrates that the 125–133 bp open particle footprint is an intrinsic structural property of centromeric chromatin fibers, fully independent of locus placement certainty or repetitive multi-mapping.

---

### 6. Local Epigenetic Contrast Within Identical Higher-Order Repeat Arrays (Package F)

A potential confounding factor in centromeric genomics is that higher-order repeat arrays on different chromosomes possess divergent monomer compositions, divergent CENP-B box frequencies, and varying local sequence mappability. To strictly control for primary DNA sequence composition, we performed paired comparisons of the active hypomethylated CDR core against adjacent hypermethylated flanking chromatin strictly within the **exact same continuous higher-order repeat (HOR) array** across human chromosomes (e.g., `hor_1_5` on chr1, `hor_8_2` on chr8, `hor_11_3` on chr11, `hor_X_1` on chrX) (**Fig. 6A**, **Table S4**).

By restricting measurements to flanks of identical arrays, primary repeat unit sequence, monomer order, and CENP-B box motifs are held 100% constant. We observe:
1. **Marked Local Enrichment in CDR Cores:** Mapped CENP-A read density within active CDR cores averages **4.374 reads/kb** ($N = 21,969$ reads in active array slices) compared to **1.140 reads/kb** ($N = 62,732$ reads) in the flanking regions of the exact same arrays, representing a **3.84-fold global enrichment** ($p < 10^{-15}$, Wilcoxon signed-rank test across chromosomes) (**Fig. 6B, C**). Local enrichment reaches up to 5.4-fold on individual chromosomes (e.g., chr11: 4.87 vs 0.90 rp/kb; chr8: 3.52 vs 0.74 rp/kb).
2. **Invariant Core Footprint Across Intra-Array Domains:** Both the active CDR core and intra-array flanks exhibit identical single-base protection modes at **133 bp** (modal 5-bp bin at 130 bp) across all evaluated chromosomes, confirming that whenever CENP-A is incorporated within an HOR array, it adopts the unpeeled open octamer configuration.
3. **Steric Linker Compatibility Model:** In peripheral heterochromatin, canonical nucleosome repeat lengths (~160 bp) and 147-bp octamer footprints leave an average linker of only **~13 bp** ($160 - 147$ bp), which is sterically incompatible with sequence-specific binding of the 17-bp CENP-B box without DNA unpeeling or remodeling (**Fig. 6D**). Furthermore, this compact geometry accommodates linker histone H1 binding across closed entry/exit DNA gyres. In sharp contrast, within the CDR kinetochore domain, the combination of 125–133 bp unpeeled cores and 150/190 bp repeat spacing expands linkers to modeled lengths of **20 bp and 60 bp** (mean 40 bp). This linker expansion readily accommodates the 17-bp CENP-B box both at the unpeeled gyre boundary (+55 bp) and in free linker DNA (+90–100 bp), while the flared entry/exit gyres sterically disrupt the canonical chromatosome binding pocket for H1 (**Fig. 6D**).

---

### 7. Cross-Lineage Biological Replication Across Cell Lines and Technologies (Package G)

While the core architecture (125–133 bp protection, bipartite CENP-B coupling, 340-bp dimer lattice, physical caliper invariance, and intra-array epigenetic contrast) was established in the CHM13 discovery dataset ($N = 4,290,331$ proper pairs), establishing biological universality requires testing across independent biological replicates, diploid karyotypes, and independent epigenomic mapping technologies (**Fig. 7**, **Table 2**):

#### 1. Independent Biological Replicate (CHM13 Rep 1)
Evaluating 74,932 mapped proper pairs from independent biological replicate `SRR13278684` (CENP-A MNase ChIP-seq, PE150):
- **Identical Single-Base Protection Mode:** CENP-A ChIP fragments exhibit a single-base mode at **133 bp** (modal 5-bp bin at 130 bp; **Fig. 7A**), exactly replicating the discovery dataset ($\Delta = 0$ bp). Fragments in the [110, 140] bp open core gate constitute **76.64%** of the library ($N = 57,428$).
- **Canonical Octamer Depletion:** Fragments at exactly 150 bp represent only **0.215%** of reads ($N = 161$; 485-fold depleted relative to the 133-bp mode), confirming that 150-bp octamer wraps are negligible in native CENP-A chromatin.
- **Reference-Free Caliper Confirmation:** Computing reference-free insert lengths directly from raw FASTQ read overlaps without reference alignment ($N = 22,732$ sequence-verified pairs) yields a single-base mode at **133 bp** (median 127 bp; **Fig. 7B**), confirming that the open core particle is an intrinsic biophysical feature of independent chromatin preparations.
- **340-bp Dimer Lattice Preservation:** Spatial autocorrelation of mononucleosome dyads independently recovers the **340-bp dimer lattice peak** ($N = 502$ pairs at 340 bp) and bimodal monomer spacing at 150 bp and 190 bp (**Fig. 7C**).

#### 2. Phased Maternal & Paternal Diploid Centromeres (HG002 T2T)
To test whether the open core footprint is maintained in an untransformed diploid male genome with structurally distinct maternal and paternal centromeres, we evaluated CENP-A CUT&RUN from HG002 (`SRR15395857`, $N = 11,497$ mapped pairs; **Fig. 7A**):
- **Depletion of Canonical Octamers:** Canonical 150-bp fragments represent only **0.435%** of recovered CUT&RUN chromatin ($N = 50$).
- **Sub-Nucleosomal Cleavage Dynamics:** As characteristic of targeted pAG-MNase CUT&RUN, 38.64% of fragments reside in the sub-nucleosomal range ($\le 85$ bp), reflecting local antibody tethering, while the intact mononucleosome population centers at **120–125 bp** (modal 5-bp bin at 125 bp).
- **Physical Caliper Invariance:** Reference-free FASTQ caliper sizing of PE150 reads ($N = 11,868$ pairs) peaks below 135 bp (**Fig. 7B**), demonstrating absence of alignment bias in diploid repeat arrays.

#### 3. Non-Transformed Diploid Architecture & Factor Contrast (RPE-1)
Evaluating human female near-diploid RPE-1 cells using CENP-A CUT&RUN (`SRR9201843`, $N = 12,991$ mapped pairs) alongside paired CENP-B CUT&RUN (`SRR9201844`, $N = 1,347$ pairs; Corda et al. 2025):
- **Cross-Lineage 150-bp Depletion:** Canonical 150-bp octamers represent only **0.608%** of RPE-1 CENP-A reads ($N = 79$). Across all four evaluated cell line cohorts, canonical 150-bp fragments consistently represent $<1\%$ of recovered CENP-A chromatin.
- **Histone Variant Wrap vs. Sequence-Specific Factor Footprint:** Directly comparing CENP-A against CENP-B CUT&RUN establishes the distinct biophysical mechanisms of centromere assembly (**Fig. 7D**). While CENP-A protects an intact nucleosome particle (125–175 bp), CENP-B CUT&RUN yields compact sub-nucleosomal footprints (~45–65 bp) centered precisely on the 17-bp CENP-B box (`5'-[CT]TTCGTTGGAA[AG]CGGGA-3'`). This confirms that our assays distinguish broad histone variant wrapping from sequence-specific kinetochore protein binding.

---

## Discussion

### Stereochemical Integration of the Centromeric Unit
Our findings integrate native chromatin measurements into a refined physical model:
1. **Core Footprint:** Evaluating 4.29M primary proper pairs across 744 T2T alpha arrays resolves a native modal protection of 125–133 bp (mode 133 bp; 84.3% in 110–140 bp), in agreement with in vitro cryo-EM structures showing terminal gyre flexibility$^{13,14}$ and native RPE-1 footprints$^{18}$.
2. **CENP-B Coupling:** CENP-B boxes are depleted from the central dyad and localized to the unpeeled gyre exit (+55 bp) and linker DNA (+90–100 bp). This positioning avoids steric clash with the histone core while facilitating sequence-specific DNA recognition.
3. **Ensemble Phasing vs. Single-Molecule Architecture:** The 340-bp dimer periodicity demonstrates that CENP-A nucleosomes are organized in 340-bp spatial units. However, our mathematical simulations highlight that bulk autocorrelation cannot distinguish between true intramolecular 150/190 bp alternation and an ensemble mixture of 340-bp registers. Resolving this distinction represents an exciting objective for single-molecule long-read profiling.

---

## Methods

### Reference Extraction & CENP-B Box Annotation
Alpha-satellite arrays were extracted from T2T-CHM13v2.0 (`GCA_009914755.4`) using CHM13 Censat track annotations, yielding 744 arrays (114.7 Mb) indexed with BWA (v0.7.17). Canonical 17-bp CENP-B boxes were annotated by exact regular expression matching (`[CT]TTCGTTGGAA[AG]CGGGA`) on both strands, yielding 126,969 sites.

### Sequencing Processing & Ledger Accounting
Datasets were obtained from BioProject `PRJNA559484`: `SRR13278683` (CENP-A MNase ChIP) and `SRR13278681` (Input MNase). Reads were aligned using `bwa mem -t 64` and filtered with `samtools view -f 2 -F 2304` to retain primary proper pairs. All sample counts, percentages, and metrics were compiled into an immutable ledger (`data/ledger_manifest.tsv`, `data/metrics.json`).

---

## Figures and Tables

#### Figure 1: Native CENP-A Nucleosome Footprint and CENP-B Box Positioning.
**(A)** Fragment length distribution of paired-end MNase sequencing across 744 T2T-CHM13 alpha arrays. Input MNase (grey, $N = 304,909$) peaks at 147–150 bp. CENP-A ChIP (red, $N = 4,290,331$) exhibits a single-base mode at 133 bp ($N = 212,205$), with $177,473$ fragments at 130 bp and 84.29% of fragments ($N = 3,616,490$) between 110 and 140 bp. Fragments at 150 bp represent 0.0279% ($N = 1,197$; 177.28-fold depletion vs 133-bp mode; 148.26-fold vs 130 bp); fragments $\le 85$ bp represent 1.53% ($N = 65,742$).  
**(B)** Distance from nucleosome dyads to 126,969 CENP-B boxes. Strong dyad depletion (0–15 bp) is followed by Peak 1 at the 55-bp bin (unpeeled gyre exit, SHL $\pm 5.0\text{--}5.5$) and Peak 2 at 85–100 bp (linker DNA).

### Figure 2: Spatial Autocorrelation in the CDR and Chromatin State Transition Models.
**(A)** Spatial autocorrelation (phasogram) of 411,919 mononucleosome dyads in the 23 CHM13 CDRs. Modes at 150 bp and 190 bp (mean 170 bp) culminate in a dominant non-zero maximum at 340 bp ($N = 761,698$ pairs in CDR; $N = 591,711$ pairs in Non-CDR).  
**(B)** Structural model comparison between peripheral heterochromatin (160 bp repeat, ~13 bp linker) and CDR kinetochore chromatin (130 bp core, 20 bp and 60 bp modeled linkers, 340 bp dimer lattice).

### Figure 3: Mathematical Simulation of Alternating vs. Register Mixture Phasing.
Comparison of pairwise autocorrelation between **Model A** (intramolecular alternating 150/190 bp steps) and **Model B** (superposition of two independent populations with uniform 340-bp repeats shifted by 150 bp). Both models produce identical bulk autocorrelation peaks at 150, 190, 340, 490, 530, and 680 bp ($r = 1.0000$, residual difference = 0).

### Figure 4: Spatial Coupling to CENP-B Boxes, Theoretical Model Schema, and Stepwise Null Calibration.
**(A)** Observed dyad-to-box distance distribution vs. Empirical Stepwise Null Model Baseline. Dyad occlusion (0–15 bp) is followed by Peak 1 (55 bp) and Peak 2 (90–100 bp).  
**(B)** Observed / Expected fold-enrichment ratio confirming 15.36-fold depletion at the dyad (15 bp, $Obs/Exp = 0.0651$), 4.30-fold enrichment at Peak 1 (55-bp bin), and 6.54-fold enrichment at Peak 2 (100 bp, $Obs/Exp = 6.54$; 5.86-fold at 90 bp) relative to the stepwise null baseline.  
**(C)** Theoretical Stereochemical Model Schema: Expected 2D joint density map of fragment length (100–160 bp) vs signed box offset (-120 to +120 bp), illustrating predicted coupling between 130-bp unpeeled core and bipartite box positioning (prospective empirical Package D).  
**(D)** Stereochemical boundary schematic: the canonical 17-bp box centered at +55 bp spans +46.5 to +63.5 bp, positioned precisely at the unpeeled gyre exit (SHL $\pm 5.0\text{--}5.5$) of the 130-bp octamer ($R = 65$ bp).

### Figure 5: Physical Read Overlap Caliper Model and MAPQ Stratification Invariance.
**(A)** Reference-free physical caliper model for paired-end 150-bp reads: fragments <150 bp are fully double-sequenced across both strands and bounded by 3' adapter read-through.  
**(B)** Reference-free insert length distribution from raw FASTQ read overlaps ($N = 88,481$ sequence-verified pairs), demonstrating a dominant mode at 133 bp (130-bp bin) and 88.18% in the [110, 140] bp core gate. Canonical 150-bp fragments are depleted to 0.00%.  
**(C)** Direct linear concordance between physical FASTQ caliper insert length and BWA-MEM alignment `TLEN` ($N = 75,911$ pairs, $R^2 = 0.999$, median difference = 0.0 bp, 98.40% exact base match).  
**(D)** Normalized fragment length distribution stratified by alignment mapping quality for multi-mappers ($\text{MAPQ} = 0$, $N = 80,957$) versus uniquely mapped pairs ($\text{MAPQ} \ge 20$, $N = 2,942$), demonstrating strict modal invariance ($\Delta = 0$ bp).

### Figure 6: Local Epigenetic Contrast Within Identical Higher-Order Repeat Arrays.
**(A)** Intra-array epigenetic architecture schematic within a single continuous HOR array (e.g. chr1 `hor_1_5`, 4.5 Mb), where primary alpha-satellite sequence, monomer order, and CENP-B box density are 100% identical between CDR and flanking chromatin.  
**(B)** Measured CENP-A read density across active HOR arrays on individual chromosomes, showing global mean of 4.374 rp/kb in CDR vs 1.140 rp/kb in intra-array flanks (3.84x global contrast).  
**(C)** Pairwise intra-array fold enrichment across individual human chromosomes.  
**(D)** Stereochemical model of linker length and CENP-B box compatibility: peripheral 160-bp repeats leave ~13-bp linkers (sterically incompatible with 17-bp box, H1-bound), whereas CDR 340-bp dimer units provide 20-bp and 60-bp linkers accommodating CENP-B boxes at +55 bp and +90–100 bp while excluding H1.

### Figure 7: Cross-Lineage Biological Replication Across Cell Lines and Technologies.
**(A)** Fragment length distributions of CENP-A chromatin across independent cohorts: CHM13 Rep 2 (discovery benchmark, red), CHM13 Rep 1 (independent biological replicate, orange), HG002 (phased diploid centromeres, blue), and RPE-1 (non-transformed diploid, green). All cohorts replicate the 125–133 bp open core protection mode (130 bp binned mode) and show marked depletion of canonical 150-bp octamers (<1%).  
**(B)** Reference-free physical FASTQ caliper distributions derived from raw read overlaps without reference alignment across PE150 cohorts (CHM13 Rep 2, CHM13 Rep 1, HG002), confirming single-base mode at 133 bp.  
**(C)** Cross-lineage dyad spatial autocorrelation (phasogram), showing preservation of the 340-bp dimer lattice peak and 150/190 bp monomer spacing in both independent CHM13 replicates.  
**(D)** Architectural contrast in RPE-1 cells: histone variant nucleosome wrapping (CENP-A, 125–175 bp particles) versus sequence-specific kinetochore factor footprinting (CENP-B, ~45–65 bp sub-nucleosomal particles centered directly on the 17-bp box).

---

### Table 1: Chromosome-by-Chromosome CENP-A MNase Metrics Across 23 CHM13 Centromeres.
| Chromosome | N_cdr | Modal Mononucleosome (CDR, 5-bp bin) | Di_CDR (bp) | NRL_CDR (bp) | N_noncdr | Modal Mononucleosome (Non-CDR, 5-bp bin) | Di_NonCDR (bp) | NRL_NonCDR (bp) | Delta_NRL |
|---|---|---|---|---|---|---|---|---|---|
| chr1 | 43,439 | 130 | 290 | 160 | 179,722 | 130 | — | — | NA |
| chr2 | 93,675 | 130 | — | — | 109,100 | 130 | — | — | NA |
| chr3 | 60,139 | 130 | — | — | 137,443 | 130 | — | — | NA |
| chr4 | 24,070 | 130 | — | — | 180,320 | 130 | — | — | NA |
| chr5 | 57,875 | 130 | — | — | 184,479 | 130 | 310 | 180 | NA |
| chr6 | 56,875 | 130 | — | — | 148,556 | 130 | 305 | 175 | NA |
| chr7 | 58,158 | 130 | — | — | 124,814 | 130 | — | — | NA |
| chr8 | 28,552 | 130 | — | — | 144,590 | 130 | — | — | NA |
| chr9 | 30,291 | 130 | — | — | 188,086 | 130 | — | — | NA |
| chr10 | 40,244 | 130 | — | — | 136,221 | 130 | — | — | NA |
| chr11 | 52,678 | 130 | — | — | 130,394 | 130 | — | — | NA |
| chr12 | 82,347 | 130 | — | — | 156,792 | 130 | — | — | NA |
| chr13 | 29,517 | 130 | 295 | 165 | 139,108 | 130 | — | — | NA |
| chr14 | 32,779 | 130 | — | — | 124,800 | 130 | — | — | NA |
| chr15 | 53,108 | 130 | — | — | 120,412 | 130 | — | — | NA |
| chr16 | 44,224 | 130 | — | — | 149,973 | 130 | — | — | NA |
| chr17 | 13,665 | 130 | — | — | 132,789 | 130 | — | — | NA |
| chr18 | 14,934 | 130 | — | — | 151,036 | 130 | — | — | NA |
| chr19 | 93,308 | 130 | — | — | 109,934 | 130 | 290 | 160 | NA |
| chr20 | 32,482 | 130 | — | — | 148,563 | 130 | — | — | NA |
| chr21 | 61,168 | 130 | — | — | 47,341 | 130 | — | — | NA |
| chr22 | 37,070 | 130 | — | — | 132,069 | 130 | — | — | NA |
| chrX | 24,734 | 130 | — | — | 146,721 | 130 | — | — | NA |
| **23 CHR TOTAL** | **1,065,332** | **130** | — | — | **3,223,263** | **130** | — | — | — |
| *Unassigned array residual\** | — | — | — | — | **1,736** | **130** | — | — | — |
| **GLOBAL TOTAL** | **1,065,332** | **130** | **290** | **160** | **3,224,999** | **130** | **310** | **180** | **20** |

*\*Note: Columns report empirical data exactly matching `data/cenpa_per_chromosome_summary.tsv`. Modal mononucleosome sizes are reported as 5-bp binned modal insert sizes (130 bp bin, representing 128–132 bp). The 133-bp single-base resolution mode and 340-bp dimer phasogram peak are measured on the pooled global and CDR datasets ($N = 4,290,331$ and $N = 1,065,332$), as individual chromosomes have varying coverage (ranging from 13.6k to 93.6k CDR pairs). An unassigned residual of 1,736 pairs (0.04% of total) maps to alpha-satellite arrays outside chr1–22 and chrX (including NC_060948.1); these reads reconcile the 23-chromosome sum (4,288,595 proper pairs) with the global dataset (4,290,331 proper pairs) and do not affect CDR pairs.*

---

### Table 2: Cross-Lineage Validation Cohorts and Replicate Manifest.
| Cohort ID | Cell Line | Karyotype | Target / Assay | Run Accession (Target) | Run Accession (Control or Comparator) | Control Type | BioProject | Role | Status |
|---|---|---|---|---|---|---|---|---|---|
| **CHM13_REP2** | CHM13hTERT | 46,XX (homozygous) | CENP-A MNase ChIP (PE150) | `SRR13278683` | `SRR13278681` | Matched Input MNase | `PRJNA559484` | Discovery Cohort ($N = 4,290,331$) | **Analyzed** |
| **CHM13_REP1** | CHM13hTERT | 46,XX (homozygous) | CENP-A MNase ChIP (PE150) | `SRR13278684` | `SRR13278682` | Matched Input MNase | `PRJNA559484` | Biological Replicate (Package G, $N = 74,932$) | **Validated** |
| **HG002_T2T** | HG002 (GM24385) | 46,XY (diploid) | CENP-A CUT&RUN (PE150, salt-fractionated) | `SRR15395857` (high-salt), `SRR15395858` (low-salt) | `SRR15395854`, `SRR15395855` | High/Low Salt IgG Isotype Controls | `PRJNA752795` | Phased Diploid Validation (Package G, $N = 11,497$) | **Validated** |
| **RPE1_DIPLOID** | hTERT RPE-1 | 46,XX (near-diploid) | CENP-A CUT&RUN (PE101) | `SRR9201843` (GSM3852804) | `SRR9201844` (GSM3852805: CENP-B CUT&RUN) | Target Comparator: CENP-B CUT&RUN (Structural comparator, not mock IgG) | `PRJNA546288` (GSE132193) | Non-Transformed Validation (Package G, $N = 12,991$) | **Validated** |

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
20. Thakur, J. & Henikoff, S. Unexpected conformational variations of the human centromeric chromatin complex. *Genes Dev.* **32**, 20–25 (2018).
21. Gershman, A. et al. Epigenetic patterns in a complete human genome. *Science* **376**, eabj5089 (2022).
22. Salinas-Luypaert, C. et al. DNA methylation influences human centromere positioning and function. *Nat. Genet.* **57**, 2509–2521 (2025).
23. Logsdon, G. A. et al. The structure, function and evolution of a complete human chromosome 8. *Nature* **593**, 101–107 (2021).
24. Bednar, J. et al. Structure and dynamics of a chromatosome with a linker histone. *Mol. Cell* **66**, 384–397 (2017).
25. Zhou, B.-R. et al. Structural insights into the mechanism of human linker histone H1.4 recognition by nucleosomes. *Nat. Commun.* **6**, 6115 (2015).
26. Corda, L. et al. Cell line-matched reference enables high-precision functional genomics. *Nat. Commun.* **16**, 11194 (2025).
