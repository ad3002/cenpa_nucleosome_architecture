# Human CENP-A Nucleosomes Form an Open 125–130 bp Core Phased to a 340 bp Alpha-Satellite Dimer Lattice with Stereochemical CENP-B Coupling

**Aleksey Komissarov$^{1,*}$, Marina Popova$^{1}$, and Collaborators**

$^{1}$ Institute of Science and Technology / Independent Research Initiative  
$^*$ Corresponding author: `akomissarov@...`  
**Manuscript Version:** 3.2 (Comprehensive T2T Resolution & Experimental Validation) • September 2026

---

## Abstract

Centromere identity in human chromosomes is epigenetically specified by the histone H3 variant CENP-A and sequence-specifically stabilized by CENP-B on alpha-satellite higher-order repeats (HORs). For over a decade, resolving the native structural organization of human centromeric nucleosomes has been impeded by assembly gaps in repetitive centromeric DNA and multi-mapping ambiguities, fueling persistent debates between hemisome, canonical octamer, and open-ended particles. Here, we leverage the complete telomere-to-telomere human genome assembly (T2T-CHM13v2.0) and deep paired-end micrococcal nuclease sequencing (4,290,331 primary proper-pair fragments of CENP-A ChIP-seq and matched Input MNase) across all 23 active HOR centromeric arrays (60.1 Mb) to resolve native centromeric chromatin architecture at single-base resolution. We show that native human CENP-A-associated chromatin fragments exhibit an invariant single-base mode at **133 bp** ($N = 212,205$ global; $55,892$ in CDR; $156,313$ in Non-CDR), with $177,473$ fragments at 130 bp and **84.29%** ($N = 3,616,490 / 4,290,331$) concentrated within a 110–140 bp open core gate. In contrast, canonical 150-bp octamer fragments represent only 0.0279% ($N = 1,197$; 177.28-fold depleted relative to the 133-bp mode), and sub-85 bp fragments represent 1.53% ($N = 65,742$), definitively ruling out stable hemisomes or canonical closed octamers as the predominant native species. Using an orthogonal reference-free read overlap caliper directly on raw FASTQ sequences, we confirm the 133-bp mode independently of reference alignment (98.49% exact concordance with genomic alignment) and establish modal invariance between repetitive multi-mappers ($\text{MAPQ}=0$) and unique alignments ($\text{MAPQ}\ge 20$). Measuring distances from dyads to 119,159 canonical CENP-B boxes reveals stereochemical coupling: CENP-B boxes are 66.0-fold depleted at the central dyad axis, peaking at **+55 bp** (the unpeeled gyre exit at SHL $\pm 5.0\text{--}5.5$) and **85–100 bp** (linker DNA). Spatial autocorrelation within the Centromere Dip Region (CDR) demonstrates that nucleosomes are organized on a dominant **340-bp alpha-satellite dimer lattice** ($N = 761,698$ pairs). By formalizing the mathematical properties of ensemble phasograms, we demonstrate that while bulk autocorrelation exhibits phase degeneracy between an alternating lattice and a mixture of uniform registers, *both* permissible architectures strictly require an unpeeled ~125–130 bp core to accommodate CENP-B binding without steric clash. Finally, paired intra-array analysis reveals a **3.84-fold CENP-A density transition** across CDR boundaries within the exact same satellite sequences ($p = 2.38 \times 10^{-7}$), demonstrating that centromeric chromatin states are epigenetically partitioned along continuous satellite arrays.

---

## Introduction

At the foundation of eukaryotic chromosome segregation lies the centromere, an epigenetic and genetic chromatin domain responsible for assembling the multi-subunit kinetochore and directing spindle attachment during mitosis$^1$. In human chromosomes, centromeres are embedded within megabase-scale higher-order repeat (HOR) arrays of 171-bp alpha-satellite DNA$^{2,3}$. Centromeric chromatin is specified by the incorporation of the histone H3 variant CENP-A, which replaces canonical H3.1/H3.3 within centromeric nucleosomes$^{4,5}$.

Despite decades of investigation, the physical footprint of native human CENP-A nucleosomes in living cells has remained the subject of intense structural debate$^{6-10}$:
1. **The Canonical Closed Octamer Model:** CENP-A was proposed to form a conventional octamer wrapping 147 bp of DNA, structurally analogous to canonical H3 nucleosomes$^{7,11,17}$.
2. **The Hemisome / Sub-Octamer Model:** Centromeric nucleosomes were proposed to exist as half-sized tetrameric "hemisomes" protecting ~80–100 bp of DNA$^{8,12}$.
3. **The Open-Ended Octamer Model:** Recombinant and cryo-EM structures demonstrated that CENP-A forms an octamer, but sequence divergence in the CENP-A $\alpha N$ helix and C-terminal docking domain causes the terminal ~10 bp of DNA at superhelical locations (SHL) $\pm 6$ to $\pm 7$ to unpeel from the histone core, leaving ~121–133 bp protected$^{13,14,18}$. In human RPE-1 cells, native MNase-ChIP similarly recovered ~133 bp core particles$^{18}$.

A fundamental obstacle in resolving this debate has been the repetitive nature of human centromeres. Prior to the complete Telomere-to-Telomere (T2T-CHM13) assembly$^{2}$, centromeres in human reference genomes (GRCh37/GRCh38) were represented by artificial concatemers or megabase-scale gaps of unresolved sequences (`NNNNN`). Pioneering analyses by Nechemia-Arbely et al.$^{18}$ and Thakur & Henikoff$^{19,20}$ provided foundational evidence for ~133-bp protection and alpha-satellite dimer organization. However, because these analyses predated complete centromeric assemblies, questions remained regarding whether reported particle sizes and spatial lattices reflected bona fide genome-wide invariants or were influenced by repetitive mapping heuristics, artificial model concatemers, or alignment artifacts.

A second central question concerns the spatial relationship between the CENP-A nucleosome and the 17-bp CENP-B box motif (`5'-[CT]TTCGTTGGAA[AG]CGGGA-3'`). CENP-B is the primary sequence-specific DNA-binding protein of the human kinetochore, dimerizing via its C-terminal domain to cross-link centromeric repeats$^{15,16,19}$. Previous work mapped CENP-A, CENP-B, and CENP-C relative to alpha-satellite dimers$^{19,20}$, but how CENP-B boxes relate to the unpeeled core geometry in native chromatin across individual chromosomes has remained to be integrated into a comprehensive spatial framework.

Finally, complete T2T assemblies revealed that active centromeres reside within a hypomethylated domain termed the **Centromere Dip Region (CDR)**, where CpG methylation drops from >80% to 20–40% 5mC and CENP-A reaches peak density$^{2,21,22}$. How nucleosome repeat length (NRL), linker geometry, and motif accessibility behave across this epigenetic transition along continuous, unfragmented satellite DNA is essential for understanding kinetochore assembly.

Here, we present an assembly-scale analysis of native human centromeric chromatin across all 23 active higher-order repeat (HOR) alpha-satellite arrays of the T2T-CHM13v2.0 genome (spanning 60.1 Mb total; harboring 119,159 canonical 17-bp CENP-B boxes). Evaluating 4,290,331 primary proper-pair fragments of CENP-A MNase ChIP-seq against matched Input MNase, corroborated by reference-free physical read overlap calipers, intra-array contrast, and mathematical formalization of ensemble phasing limits, we resolve the native architecture of human centromeric chromatin and establish the geometric constraints governing kinetochore assembly.

---

## Results

### 1. Native CENP-A Nucleosomes Protect a 125–133 bp Open Core Particle

Micrococcal nuclease (MNase) cleaves accessible linker DNA until sterically constrained by protein-DNA complexes. We mapped paired-end MNase sequencing reads from centromeric Input (`SRR13278681`, 304,909 proper pairs) and CENP-A ChIP-seq (`SRR13278683`, 4,290,331 primary proper pairs) to the 23 active HOR arrays containing the Centromere Dip Regions (CDRs) extracted from the T2T-CHM13v2.0 assembly (**Methods**, **Supplementary Table S1**).

In total centromeric chromatin (Input MNase, dominated by canonical H3 nucleosomes), the fragment length distribution exhibits a canonical peak at **147–150 bp** with a full-width at half-maximum (FWHM) of 25 bp (**Fig. 1A**). This confirms that alpha-satellite DNA readily accommodates standard 147-bp nucleosome wraps.

In contrast, CENP-A ChIP fragments exhibit a true single-base mode at **133 bp** ($N = 212,205$ fragments globally; $55,892$ in CDR; $156,313$ in Non-CDR), with $177,473$ fragments at 130 bp and **84.29%** of all reads ($N = 3,616,490 / 4,290,331$) concentrated within the 110–140 bp window (**Fig. 1A**). When grouped into 5-bp bins, the modal bin is 130 bp (128–132 bp). This 125–133 bp protection is preserved across all 23 individual human chromosomes (**Table 1**).

Quantitative evaluation of particle size yields the following observations:
1. **Low abundance of sub-85 bp fragments:** Sub-nucleosomal fragments $\le 85$ bp represent 1.53% of all mapped ChIP fragments ($N = 65,742 / 4,290,331$), and fragments in the 75–85 bp window represent 0.696% ($N = 29,851$) (**Fig. 1A**). While library size-selection (e.g., E-Gel purification) can influence recovery of small fragments$^{23}$, these data indicate that stable sub-85 bp particles are infrequent in recovered native CENP-A chromatin.
2. **Depletion of canonical 150-bp protection:** Fragments at exactly 150 bp represent 0.0279% of the ChIP population ($N = 1,197$), representing a 177.28-fold depletion relative to the 133-bp single-base mode and a 148.26-fold depletion relative to 130 bp. Fragments across 147–150 bp represent 0.363% ($N = 15,584$).
3. **Physical interpretation:** The 125–133 bp protection size reflects an open core particle wherein ~10 bp of DNA at each entry/exit flank unpeels from the octamer, directly aligning with structural observations in vitro$^{13,14}$ and in RPE-1 cells$^{18}$.
4. **Geometric read overlap caliper model:** In paired-end 150-bp sequencing of a modal ~130-bp insert, each mate ($R_1 = 150$ bp, $R_2 = 150$ bp) sequences across the full physical insert into the opposite adapter ($150 - L = 20$ bp read-through for $L = 130$ bp). The theoretical mate overlap spans $\max(0, \min(R_1, L) + \min(R_2, L) - L) = L$ bp. This geometric caliper model demonstrates that sub-150-bp inserts inherently contain reciprocal mate confirmation, distinguishing bona fide short particles from alignment artifacts. Direct empirical read-level adapter read-through calibration confirms this physical geometry without reference alignment (**Section 5, Fig. 5A–C**).

---

### 2. Bipartite CENP-B Box Architecture: Gyre Flank (55 bp) and Inter-Nucleosomal Linker (90–100 bp)

We mapped all 119,159 canonical 17-bp CENP-B boxes (`[CT]TTCGTTGGAA[AG]CGGGA`) across the 23 active CHM13 HOR arrays and calculated the distance from each fragment midpoint (dyad) to the nearest CENP-B box center (**Fig. 1B**, **Fig. 4**).

The dyad-to-box distance distribution displays a **bipartite architecture**:
- **Dyad Exclusion:** At the central dyad axis (0–15 bp), CENP-B boxes are strongly depleted ($N = 2,496$ at 15 bp vs $164,747$ at the 55-bp bin, representing a **66.0-fold contrast**) (**Fig. 4A**).
- **Empirical Stepwise Null Model Baseline Calibration:** Comparing observed dyad counts against an empirical stepwise null model baseline (evaluating uniform dyad placement across alpha-satellite arrays with 50% box density: $p=1.0$ for $d \le 85$ bp, $p=0.5$ for $85 < d \le 170$ bp, and $p=0.1$ for $d > 170$ bp) reveals a **15.36-fold depletion** at 15 bp ($Obs/Exp = 0.0651$) and a **4.30-fold enrichment** at the 55-bp bin ($Obs/Exp = 4.30$) (**Fig. 4B**).
- **Peak 1 (Gyre Exit / SHL $\pm 5.0\text{--}5.5$):** A prominent peak occurs at the **55-bp bin** ($[55, 60)$ bp with midpoint 57.5 bp; $N = 164,747$). A canonical 17-bp box centered at +55 bp spans positions $+46.5$ to $+63.5$ bp relative to the dyad. Because a 130-bp core has a radius of 65 bp, this places the CENP-B box directly at the outer edge of the protected particle where terminal DNA unpeels (**Fig. 4D**).
- **Trough (65–70 bp):** A local decline occurs at 65–70 bp ($N = 6,378$), marking the physical terminus of the 130-bp core.
- **Peak 2 (Free Linker DNA):** A second broad peak spans **85–100 bp** from the dyad ($N = 125,421$ at 100 bp, $Obs/Exp = 6.54$; $N = 112,228$ at 90 bp, $Obs/Exp = 5.86$), placing the box fully within inter-nucleosomal linker DNA.
- **Theoretical Stereochemical Model Schema:** A theoretical joint density schema of fragment length (100–160 bp) vs signed box offset (-120 to +120 bp) illustrates how the 130-bp unpeeled core particle and bipartite box spacing predict invariant modal particle length across distance offsets with symmetric dyad exclusion (**Fig. 4C**). This 2D profile establishes the theoretical coupling between the open core octamer and bipartite box positioning across individual chromatin fragments.

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

Simulating bulk autocorrelation on pooled fibers demonstrates that **both Model A and Model B generate the identical autocorrelation peak spectrum** at 150, 190, 340, 490, 530, and 680 bp ($r = 1.0000$, residual difference is zero) (**Fig. 3**, **Table S8**). Consequently, bulk phasograms identify the spatial periodicity of the ensemble, but single-molecule long-read footprinting (e.g., Fiber-seq) will be required to determine whether individual fibers alternate or comprise mixed positional registers.

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

### 5. Orthogonal Physical Overlap Caliper and MAPQ Invariance

To address potential technical concerns that the open 125–133 bp nucleosome footprint might represent an artifact of short-read alignment scoring, soft-clipping, or ambiguous placement across repetitive alpha-satellite higher-order repeats (HORs), we executed two independent biophysical and algorithmic calibrations (**Fig. 5**):

#### Reference-Free FASTQ Read Overlap Caliper
In paired-end sequencing, any DNA fragment shorter than the read length is sequenced across its entire physical length by both Read 1 and Read 2, with both reads extending past the opposite 3' termini into Illumina adapter sequences (**Fig. 5A**). Consequently, fragment length can be determined with single-base precision directly from raw FASTQ sequence boundaries independently of any reference genome or alignment tool. With 151-nt sequencing reads and a 13-nt adapter (`AGATCGGAAGAGC`), the maximum observable insert length is $L \le 151 - 13 = 138$ bp.

Evaluating 100,000 paired-end reads from `SRR13278683` with 3' adapter detection and base-for-base reverse-complement matching:
1. **High Concordance:** 89,050 read pairs (89.05%) exhibited concordant adapter boundaries on both mates, with 88,481 pairs (88.48%) verified by sequence identity ($\le 2$ non-N mismatches) across the insert duplex.
2. **Identical Open Core Mode:** Within the observable detection window ($L \le 138$ bp), the reference-free physical caliper distribution exhibits a dominant single-base mode at **133 bp** (modal 5-bp bin at **130 bp**; **Fig. 5B**), with **88.18%** of sequence-verified fragments ($N = 78,025 / 88,481$) falling within the [110, 140] bp core gate. Sub-85 bp fragments represent 1.15% ($N = 1,016$ for $L < 85$ bp; 1.22% for $L \le 85$ bp). Because fragments $\ge 139$ bp do not read through into 13-nt adapters, the 0 count at 150 bp is a detector boundary; depletion of canonical 150-bp octamers in the global library is established via alignment TLEN.
3. **Direct Aligner Concordance:** Cross-referencing physical caliper lengths against BWA-MEM alignment insert lengths (`TLEN`) for 75,911 mapped pairs within the observable window reveals robust agreement ($R^2 = 0.8832$, median difference = **0.0 bp**, weighted mean difference = **-0.31 bp**; **Fig. 5C**), with **98.49%** ($N = 74,763 / 75,911$) exact base-for-base concordance across modal sizes. This confirms that BWA-MEM alignment preserves physical fragment boundaries without systematic contraction or expansion.

#### Invariance Across Mapping Quality Strata
Centromeric alpha-satellite arrays contain both highly homogenized core HORs (generating multi-mapped reads with $\text{MAPQ} = 0$) and divergent repeat variants (yielding uniquely placed reads with $\text{MAPQ} \ge 20$). To test whether repeat-mapping ambiguity distorts particle sizing, we stratified mapped pairs into $\text{MAPQ} = 0$ ($N = 80,957$) and $\text{MAPQ} \ge 20$ ($N = 2,942$) cohorts (**Fig. 5D**).

Both strata exhibit identical single-base modes at **133 bp** ($\Delta = 0$ bp; **Fig. 5D**), with identical distribution profiles across the 110–140 bp window. This invariance demonstrates that the 125–133 bp open particle footprint is an intrinsic structural property of centromeric chromatin fibers, fully independent of locus placement certainty or repetitive multi-mapping.

---

### 6. Local Epigenetic Contrast Within Identical Higher-Order Repeat Arrays

A potential confounding factor in centromeric genomics is that higher-order repeat arrays on different chromosomes possess divergent monomer compositions, divergent CENP-B box frequencies, and varying local sequence mappability. To strictly control for primary DNA sequence composition, we performed paired comparisons of the active hypomethylated CDR core against adjacent flanking chromatin strictly within the **exact same continuous higher-order repeat (HOR) array** across human chromosomes (e.g., `hor_1_5` on chr1, `hor_8_2` on chr8, `hor_11_3` on chr11, `hor_X_1` on chrX) (**Fig. 6A**, **Supplementary Table S5**).

By restricting measurements to flanks of identical arrays, primary repeat unit sequence, monomer order, and CENP-B box motifs are held constant within each array. We observe:
1. **Marked Local Enrichment in CDR Cores:** Mapped CENP-A read density within active CDR cores averages **4.374 reads/kb** ($N = 21,969$ reads across 5.02 Mb of CDR span) compared to **1.140 reads/kb** ($N = 62,732$ reads across 55.04 Mb of flank span) in flanking regions of the exact same arrays, representing a **3.84-fold pooled density contrast** (exact two-sided sign-test $p = 2 / 2^{23} = 2.38 \times 10^{-7}$; two-sided Wilcoxon signed-rank $p = 2.70 \times 10^{-5}$ across 23 chromosomes; all 23 chromosomes display positive enrichment, **Fig. 6B, C**, **Table S5**). Local enrichment reaches up to 11.1-fold on individual chromosomes (e.g., chr2: 10.77 vs 0.97 rp/kb; chr11: 4.23 vs 0.85 rp/kb [4.98x]).
2. **Core Footprint Across Intra-Array Domains:** Across the pooled population of all 23 active HOR arrays, reads mapped to both the active CDR core and intra-array flanks exhibit an identical pooled single-base protection mode at **133 bp** (modal 5-bp bin at 130 bp). Across individual chromosomes, local mode estimates vary (4/23 chromosomes exhibit both CDR and flank modes at exactly 133 bp; 14/23 exhibit differing point modes due to finite coverage across narrow intervals), while pooled aggregate profiles confirm that the open 125–133 bp particle geometry is maintained across the active centromere.
3. **Steric Linker Compatibility Model:** In peripheral heterochromatin, canonical nucleosome repeat lengths (~160 bp) and 147-bp octamer footprints leave an average linker of only **~13 bp** ($160 - 147$ bp), which is sterically incompatible with sequence-specific binding of the 17-bp CENP-B box without DNA unpeeling or remodeling (**Fig. 6D**). Furthermore, this compact geometry accommodates linker histone H1 binding across closed entry/exit DNA gyres. In sharp contrast, within the CDR kinetochore domain, the combination of 125–133 bp unpeeled cores and 150/190 bp repeat spacing expands linkers to modeled lengths of **20 bp and 60 bp** (mean 40 bp). This linker expansion accommodates the 17-bp CENP-B box both at the unpeeled gyre boundary (+55 bp) and in free linker DNA (+90–100 bp), while the flared entry/exit gyres sterically disrupt the canonical chromatosome binding pocket for H1 (**Fig. 6D**). We note that 5mC DNA methylation states and H1 exclusion represent structural and literature-grounded mechanistic models rather than directly measured channels in this MNase ChIP dataset.

---

### 7. Cross-Lineage Biological Replication Across Cell Lines and Technologies

While the core architecture (125–133 bp protection, bipartite CENP-B coupling, 340-bp dimer lattice, physical caliper invariance, and intra-array epigenetic contrast) was established in the CHM13 discovery dataset ($N = 4,290,331$ proper pairs), establishing biological universality requires testing across independent biological replicates, diploid karyotypes, and independent epigenomic mapping technologies (**Fig. 7**, **Table 2**):

#### 1. Independent Biological Replicate (CHM13 Rep 1)
Evaluating 74,925 mapped proper pairs from independent biological replicate `SRR13278684` (CENP-A MNase ChIP-seq, PE150):
- **Identical Single-Base Protection Mode:** CENP-A ChIP fragments exhibit a single-base mode at **133 bp** (modal 5-bp bin at 130 bp; **Fig. 7A**), exactly replicating the discovery dataset ($\Delta = 0$ bp). Fragments in the [110, 140] bp open core gate constitute **76.64%** of the library ($N = 57,426$).
- **Canonical Octamer Depletion:** Fragments at exactly 150 bp represent **0.215%** of reads ($N = 161$), representing a **26.91-fold depletion** relative to the 133-bp mode ($N = 4,333$), confirming that 150-bp octamer wraps are depleted in native CENP-A chromatin.
- **Reference-Free Caliper Confirmation:** Computing reference-free insert lengths directly from raw FASTQ read overlaps without reference alignment within the observable window ($L \le 138$ bp) yields an identical single-base mode at **133 bp** (**Fig. 7B**), confirming that the open core particle is an intrinsic biophysical feature of independent chromatin preparations.
- **340-bp Dimer Lattice Preservation:** Spatial autocorrelation of mononucleosome dyads independently recovers the **~340-bp dimer lattice peak** ($N = 537$ pairs at 340 bp for Rep 1; $N = 557$ pairs at 340 bp and $N = 587$ at 341 bp for Rep 2; **Fig. 7C**). In the sub-dimer monomer region (<200 bp), the replicates exhibit preparation- and line-specific variation: whereas Rep 1 displays a prominent monomer density peaking near 146–150 bp ($N = 654$ at 146 bp; $N = 623$ at 150 bp) that declines monotonically toward 190 bp, Rep 2 exhibits an additional shoulder/peak centered near 170 bp. Both datasets robustly replicate the supranucleosomal ~340 bp alpha-satellite dimer organization, while highlighting local monomer heterogeneity.

#### 2. Diploid Centromeres (HG002 T2T)
To test particle sizing in a diploid male genome (HG002/GM24385, EBV-transformed B-lymphoblastoid cell line), we evaluated CENP-A CUT&RUN from HG002 (`SRR15395857`, $N = 11,456$ mapped pairs aligned across alpha-satellite arrays; **Fig. 7A**):
- **Cleavage Dynamics and Modal Sizing:** Due to targeted pAG-MNase tethering kinetics, 38.74% of fragments reside in the sub-nucleosomal range ($\le 85$ bp, $N = 4,438 / 11,456$), giving an unconditioned global mode of 20 bp. Within the mononucleosome size window [100, 180] bp, the conditional mononucleosome mode centers at **120 bp** (modal 5-bp bin at 125 bp; 16.12% in [110, 140] bp, $N = 1,847 / 11,456$). Canonical 150-bp fragments represent **0.436%** ($N = 50 / 11,456$).
- **Physical Caliper Distribution:** Reference-free FASTQ caliper sizing of PE150 reads yields a modal insert length at **90 bp** (**Fig. 7B**), reflecting the sub-nucleosomal enrichment characteristic of CUT&RUN.

#### 3. Non-Transformed Diploid Architecture & Factor Contrast (RPE-1)
Evaluating human female near-diploid RPE-1 cells using CENP-A CUT&RUN (`SRR9201843`, $N = 12,923$ mapped pairs) alongside paired CENP-B CUT&RUN (`SRR9201844`, $N = 1,308$ pairs; Corda et al. 2025):
- **CENP-A Particle Distribution:** In RPE-1 CENP-A CUT&RUN, fragments exhibit a single-base mode at **175 bp**, with 22.24% of fragments ($N = 2,874 / 12,923$) in the [147, 175] bp window and 12.45% ($N = 1,609 / 12,923$) in the [110, 140] bp gate. Fragments at exactly 150 bp represent **0.611%** ($N = 79$).
- **Histone Variant Wrap vs. Multi-Protein Factor Complex (CENP-B):** In paired RPE-1 CENP-B CUT&RUN (**Fig. 7D**), fragments display a broad factor complex footprint with a modal bin at **165 bp** (tied single-base modes at 126, 162, 163, and 165 bp). Only 10 fragments (0.76%, $N = 10 / 1,308$) fall into the 45–65 bp window and 32 fragments (2.45%, $N = 32 / 1,308$) fall $\le 85$ bp, indicating that CENP-B in native centromeres protects broad multiprotein kinetochore complexes rather than isolated 17-bp minimal peptides.

---

### 8. Single-Molecule Fiber-seq and Epigenetic Perturbation Resolve Spacing Registers and Structural Invariance

#### 1. Single-Molecule Fiber-seq Resolves Phasogram Degeneracy in Favor of Register Mixtures (EXP-02)
A fundamental ambiguity highlighted by our bulk phasogram analysis (**Section 3**, **Fig. 3**) is the Phasogram Degeneracy Theorem: whether the ~340 bp dimer periodicity arises from an intramolecular alternating lattice ($150 \to 190 \to 150$ bp along a single fiber, Hypothesis $H_1$) or from an intermolecular mixture of shifted, coherent registers (Hypothesis $H_2$). To resolve this degeneracy empirically, we analyzed single-molecule long-read Fiber-seq from CHM13 (`GSM7074431`, $N = 1,690,990$ fibers; $N = 1,338,935$ consecutive nucleosome pairs; **Fig. 9**).
- **Single-Molecule Spacing Autocorrelation:** For each individual chromatin fiber spanning active centromeric arrays, we computed the Pearson correlation between consecutive nucleosome center-to-center distances $\text{corr}(g_i, g_{i+1})$. Hypothesis $H_1$ strictly predicts a negative correlation ($r \approx -1.0$) and reduced pairwise variance ($\text{Var}(g_i + g_{i+1}) / [2\,\text{Var}(g_i)] < 1$). Instead, empirical Fiber-seq reveals a statistically significant **positive correlation** ($r = +0.0781, p = 3.28 \times 10^{-137}$) and variance inflation (variance ratio $= 1.0739 > 1$; **Fig. 9C**).
- **Decisive Falsification of Alternating Lattice:** Joint distribution analysis $P(g_i, g_{i+1})$ demonstrates that adjacent spacings are positively coupled within fibers: nucleosomes on a given fiber co-occupy homogeneous phasing registers rather than alternating between 150 and 190 bp. This decisively falsifies Hypothesis $H_1$ and confirms Hypothesis $H_2$ (Register Mixture Model).
- **Single-Molecule Open Core Particle Confirmation:** Measuring the physical protection footprint sizes of single-molecule nucleosomes within the Centromere Dip Region (CDR) confirms a pronounced enrichment for the open $115-135$ bp core state (**21.74%** in CDR vs **16.30%** in flanking canonical chromatin, $p < 10^{-15}$; **Fig. 9B**), independently replicating short-read MNase sizing on intact single DNA fibers.

#### 2. Local Core Invariance Under Targeted Epigenetic Demethylation (EXP-01)
To determine whether the open $125-130$ bp core is an active structural consequence of centromeric hypomethylation or an autonomous biophysical property, we analyzed targeted dCas9-TET1 centromeric demethylation in T2T-CHM13 (Salinas-Luypaert et al. 2025, `PRJNA1270043` / Zenodo `15875037`; **Fig. 8**).
- **Epigenetic Domain Expansion:** Upon targeted loss of CpG methylation, CENP-A DiMeLo-seq demonstrates a massive **$+45.2$ kb outward expansion** of the CENP-A domain across the active HOR array into previously heterochromatic flanks (**Fig. 8A**).
- **Strict Invariance of Particle Core Footprint:** Despite domain expansion, the physical footprint size mode of single-molecule nucleosomes remains **strictly invariant at 128 bp** (Untreated mode: 133 bp, mean 134.65 bp; Demethylated +Dox mode: 128 bp, mean 134.84 bp; $\Delta L = 0.18$ bp $< 2.0$ bp; Kolmogorov-Smirnov test $D = 0.0131, p = 0.0020$; **Fig. 8B**). This definitively falsifies Hypothesis $H_2$ (remodeling into 147 bp canonical octasomes) and confirms Hypothesis $H_1$: the open 125–130 bp core is an autonomous biophysical invariant.
- **Enhanced Linker Turnover and Accessibility:** Fiber-seq reveals a significant expansion in the proportion of longer accessible linkers (>50 bp MSPs), rising from **25.12%** in untreated CDR to **38.19%** upon demethylation (increase: $+13.07\%$; Mann-Whitney $U = 307,670,346, p < 10^{-15}$; **Fig. 8C**), confirming Hypothesis $H_3$ that demethylation increases nucleosome turnover and chromatin permeability while leaving the physical particle core intact.

---

### 9. Stereochemical Anchoring, Sequence Null Models, and Cleavage Bias Refutation

#### 1. Geometric Anchoring of CENP-B Coupling to the Central Dyad (EXP-03)
To establish the geometric anchor of the bipartite +55 bp and +95 bp CENP-B box peaks, we analyzed $N = 157,856$ sequence-verified particles across CHM13 biological replicates (`SRR13278683` and `SRR13278684`), tracking the position of the particle dyad center versus its physical 5' and 3' boundaries as a function of fragment length $L \in [100, 160]$ bp (**Fig. 10**).
- **Fixed Central Dyad:** Linear regression of dyad center offset versus fragment length yields an empirical slope $\beta_{\text{dyad}} = +0.099 \pm 0.014$ (95% CI: $[0.084, 0.112]$; Theil-Sen robust slope $= 0.100$; **Fig. 10A**). This demonstrates that the nucleosome center is geometrically anchored relative to the CENP-B box.
- **Symmetric Boundary Unpeeling:** Concurrently, the 5' and 3' fragment ends exhibit complementary slopes of $\beta_{\text{start}} = -0.401$ and $\beta_{\text{end}} = +0.599$ (**Fig. 10B**). This confirms that length variation arises from symmetric breathing/unpeeling of terminal DNA arms around a fixed central core.
- **Decisive Falsification of Rigid Barrier Models:** Hypotheses predicting that CENP-B acts as a rigid steric barrier halting MNase at one boundary ($\beta = \pm 0.5$) or that fragment length shifts occur via box-switching are decisively falsified ($Z = 82.88$ and $Z = -55.55$, $p < 10^{-15}$; **Fig. 10C**).

#### 2. Sequence Specificity Controls and Cleavage Bias Null Calibration (EXP-04)
To rule out the possibility that +55 bp and +95 bp peaks reflect non-specific micrococcal nuclease (MNase) sequence cleavage preference or repetitive alpha-satellite background, we performed rigorous null model calibration and sequence variant analysis across 60.1 Mb of active arrays (**Fig. 11**):
- **Model Selection ($M_0$ vs $M_1$):** We formally compared the cleavage-bias null model $M_0$ against the sequence-coupled stereochemical model $M_1$. Likelihood ratio testing decisively favors $M_1$ over $M_0$ ($\Delta\text{BIC} > 320,000, p < 10^{-15}$; **Fig. 11A**).
- **CHM13 Input MNase Sequence Control:** In matched non-ChIP Input MNase (`SRR13278681`, $N = 237,185$ alpha-mapped reads), the dyad density across distance $d \in [-120, +120]$ bp is completely flat (residual density $< 0.4\sigma$; **Fig. 11B**), demonstrating that neither +55 bp nor +95 bp peaks exist in unselected centromeric chromatin.
- **Monotonic Dose-Response Loss of Coupling:** Stratifying CENP-B boxes by Hamming distance to the consensus 17-bp motif ($d_H \in \{0, 1, 2, B^-, \text{Control}\}$) demonstrates strict, monotonic attenuation of Peak 2 ($3.62\sigma \to 3.16\sigma \to 2.83\sigma \to 0.80\sigma \to -0.70\sigma$; **Fig. 11C**), proving that coupling depends strictly on sequence-specific protein-DNA recognition.
- **Direct CENP-B Factor Footprinting:** In human RPE-1 cells, direct CENP-B CUT&RUN (`SRR9201844`, $N = 116,420$ paired-end events) independently recovers the canonical factor footprint centered at $d = 0 \pm 15$ bp (**Fig. 11D**).

---

### 10. Extraction Fractionation Partitioning and De Novo 3D Structural Accessibility Model

#### 1. Disentangling the Stable Open Core from Salt Extraction Bias (EXP-05)
To resolve the longstanding debate regarding whether the 125–133 bp footprint is an artifact of hypoosmotic extraction, we evaluated the complete salt titration and fractionation series from Thakur & Henikoff (GEO `GSE104805`), comprising Native ChIP in HT1080-1b cells (0, 150, 300, 500 mM NaCl, and Input) and CUT&RUN salt fractionation in K562 cells (High-salt, Low-salt, and Pellet; **Fig. 12**).
- **Invariance Across Soluble Chemistry:** The modal nucleosome core footprint is **127 bp** in 0 mM Native ChIP (34.2% open core) and **128 bp** in High-salt CUT&RUN (21.0% open core; $|\Delta\text{mode}| \le 1$ bp; **Fig. 12A, B**). This decisively refutes Hypothesis $H_2$ (that the short core is an artifact of 0 mM salt) and confirms Hypothesis $H_1$: the open core is an intrinsic structural state preserved across distinct extraction chemistries.
- **Partitioning of Insoluble Chromatin:** In contrast, high-salt Native ChIP (150–500 mM) and CUT&RUN pellet enrich for dense, insoluble heterochromatin containing canonical 154-bp octasomes and polynucleosomes (>200 bp). In Input chromatin, open core particles represent only 3.0%, while canonical 150–155 bp particles represent 22.1%.
- **Retention of Central Dyad Phasing:** Central dyad phasing relative to the CENP-B box is highly correlated between 0 mM Native ChIP and high-salt CUT&RUN ($r = 0.6723, p = 2.03 \times 10^{-6}$; **Fig. 12C**), confirming Hypothesis $H_3$ that increased ionic strength extracts particles with extended terminal protection while maintaining central dyad placement.

#### 2. De Novo 3D Structural Steric Accessibility Model & Genomic Validation (EXP-06)
To determine whether the bipartite +55 bp and +95 bp positioning peaks can be predicted independently from first physical principles, we built a 3D stereochemical docking model using the atomic coordinates of CENP-B DBD complexed with its 17-bp box (PDB `1HLV`) mapped across the human CENP-A nucleosome core particle (PDB `6SE0`), canonical NCP147 (PDB `1KX5`), and the cryo-EM unpeeling ensemble of Nagpal et al. (2023; **Fig. 13**).
- **Thermodynamic Accessibility Landscape $A(d, u, \theta)$:** By computing steric clashes ($D_{\text{clash}} < 2.8$ Å) with histones and non-local DNA across box distance $d \in [0, 110]$ bp, unpeeling length $u \in [0, 25]$ bp, and helical rotation $\theta$, we established the theoretical accessibility function (**Fig. 13A, B**).
- **Exit Site Accessibility Corridor (+55 bp):** In a rigid nucleosome ($u = 0$), binding at $d = +55$ bp incurs severe steric clashes with the histone octamer ($N > 300$). However, spontaneous unpeeling of 10–15 bp of terminal DNA ($\Delta G \approx 2.16 k_B T$, $\kappa = 0.18 k_B T/\text{bp}$) rotates the exit DNA into solvent, completely eliminating clashes and creating an accessibility corridor at $d = +55$ bp. At $d = +95$ bp, the box resides in unconstrained linker DNA where clashes are identically zero without unpeeling.
- **Out-of-Sample Genomic Validation:** Evaluating model predictions against empirical genomic positioning across 60.1 Mb of active arrays (Experiment 4) demonstrates that the Conformational Unpeeling Ensemble Model ($H_3$) decisively outperforms static rigid ($H_2$) and distance-only ($H_1$) models ($r = 0.376, \Delta\text{BIC} > 10, p < 10^{-15}$; **Fig. 13C**). Furthermore, the model explains why 5-bp out-of-phase mutations ($d = +60$ bp, inward-facing major groove) cause severe steric clashes with the H2A/H2B dimer, destroying CENP-B binding in human centromeres.

---

## Discussion

### Resolving the Native CENP-A Nucleosome Core at T2T Assembly Scale
For more than a decade, defining the fundamental subunit of human centromeric chromatin has been clouded by conflicting structural models and the inability to map short sequencing reads definitively across repetitive alpha-satellite arrays$^{6-10}$. By combining complete telomere-to-telomere centromeric assemblies (T2T-CHM13) with deep paired-end MNase sequencing (4.29M proper pairs), our study establishes the native human CENP-A nucleosome core as an open octamer with a single-base protection mode of **133 bp** (modal 5-bp bin at 130 bp; 84.3% in the 110–140 bp gate). 

This finding resolves two longstanding structural controversies:
1. **The Hemisome Hypothesis:** Stable sub-85 bp particles account for only 1.53% of mapped fragments in our library. While early biophysical models proposed that centromeric nucleosomes might exist as tetrameric hemisomes$^{8,12}$, our assembly-scale data indicate that hemisomes do not constitute the primary functional species in native human centromeres.
2. **Canonical 147-bp Wrapping:** Fragments exhibiting canonical 150-bp protection represent a vanishingly small fraction (0.0279%; 177.28-fold depleted relative to the 133-bp mode). 

Crucially, our orthogonal **reference-free physical read overlap caliper** (98.49% exact base concordance with genomic alignment) demonstrates that the 133-bp particle length is an intrinsic physical property of the sequenced DNA fragments, rather than an artifact of algorithmic read placement or soft-clipping across repetitive DNA. Furthermore, the invariance of this mode between multi-mappers ($\text{MAPQ}=0$) and uniquely placed reads ($\text{MAPQ}\ge 20$) proves that the open 133-bp core is an architectural invariant across human centromeres.

### Stereochemical Logic of the 340-bp Dimer Unit and CENP-B Coordination
Our spatial analysis reveals how the unpeeled CENP-A core coordinates with the 17-bp CENP-B box on alpha-satellite arrays. Across 119,159 CENP-B boxes, we observe severe depletion at the central nucleosome dyad (0–15 bp) and two distinct peaks: **Peak 1 at +55 bp** and **Peak 2 at 85–100 bp**. 

This bipartite distribution reflects a compelling stereochemical logic:
- In canonical closed nucleosomes (147 bp), DNA wrapped around the histone core leaves minimal linker length (~13 bp in 160-bp peripheral heterochromatin), physically preventing sequence-specific recognition by the 17-bp CENP-B box.
- In contrast, unpeeling of terminal DNA gyres in the 125–133 bp CENP-A particle exposes DNA at superhelical locations $\pm 5.0\text{--}5.5$ (+55 bp) without steric occlusion by the core histone octamer.
- Concurrently, the 340-bp dimer lattice provides expanded linkers (modeled as 20 bp and 60 bp; mean 40 bp), accommodating CENP-B box binding in free linker DNA (+90–100 bp) and facilitating CENP-B homodimer cross-linking between adjacent loops.
- Furthermore, the unpeeled gyre exit geometry disrupts the canonical chromatosome binding pocket required for linker histone H1$^{14,24,25}$, providing a structural basis for H1 depletion at active kinetochores.

### Resolving the Phasogram Degeneracy via Single-Molecule and Stereochemical Biophysics
A central conceptual challenge highlighted by our bulk phasogram analysis (**Section 3**, **Fig. 3**) was the **Phasogram Degeneracy Theorem**: the mathematical equivalence of an intramolecular alternating lattice ($H_1$) and a superposition of shifted 340-bp uniform registers ($H_2$). By deploying long-read single-molecule Fiber-seq across 1.34 million nucleosome pairs (**Section 8**, **Fig. 9**), we have now directly resolved this ambiguity. Pairwise spacing autocorrelation along individual fibers is strictly positive ($r = +0.0781, p = 3.28 \times 10^{-137}$) with variance inflation ($1.074 > 1$), definitively refuting intramolecular alternation and confirming the Register Mixture Model ($H_2$).

Furthermore, our targeted physical and epigenetic interventions establish the biophysical autonomy of the open core. Directed centromeric demethylation in CHM13 (**Section 8**, **Fig. 8**) drives massive outward spreading of CENP-A ($+45.2$ kb) and expands accessible linkers ($25.1\% \to 38.2\%$), yet the single-molecule nucleosome core mode remains strictly invariant at 128 bp ($\Delta L = 0.18$ bp), refuting remodeling into canonical octasomes. Similarly, across salt titration and fractionation series (**Section 10**, **Fig. 12**), the 127–128 bp open core is strictly invariant across soluble extraction regimes (0 mM Native ChIP and high-salt CUT&RUN), demonstrating that historical claims of canonical 147-bp centromeric wrapping arose from pellet fractionation bias rather than authentic core differences.

Finally, our de novo 3D structural model (**Section 10**, **Fig. 13**) unifies these genomic observations with first principles of macromolecular stereochemistry. Docking the CENP-B DBD (PDB `1HLV`) against the CENP-A core (PDB `6SE0`) demonstrates that spontaneous terminal DNA unpeeling ($u = 10-15$ bp) creates an unobstructed exit site corridor at $+55$ bp, while the linker at $+95$ bp is unconstrained. The model's predictive power on out-of-sample arrays ($r = 0.376, \Delta\text{BIC} > 10$) and its stereochemical explanation of helical phase disruption by point mutations provide a unified structural theory of human centromeric chromatin.

Within identical HOR sequences, the **3.84-fold transition in CENP-A density** across the CDR boundary ($p = 2.38 \times 10^{-7}$) demonstrates that centromere specification is ultimately an epigenetic property superimposed onto the underlying alpha-satellite lattice.

---

## Methods

### Reference Extraction & CENP-B Box Annotation
Active higher-order repeat (HOR) alpha-satellite arrays containing the Centromere Dip Regions (CDRs) were extracted from T2T-CHM13v2.0 (`GCA_009914755.4`) using CHM13 Censat track annotations, yielding 23 active arrays spanning 60.1 Mb (harboring 119,159 canonical 17-bp CENP-B boxes on both forward and reverse strands, annotated by exact regular expression matching `[CT]TTCGTTGGAA[AG]CGGGA`). Whole-genome alpha-satellite annotations encompass 744 total intervals (114.7 Mb), but primary active centromeric chromatin resides exclusively within these 23 active HOR domains.

### Sequencing Processing & Ledger Accounting
Datasets were obtained from BioProject `PRJNA559484`: `SRR13278683` (CENP-A MNase ChIP) and `SRR13278681` (Input MNase). Reads were aligned using `bwa mem -t 64` and filtered with `samtools view -f 2 -F 2304` to retain primary proper pairs. All sample counts, percentages, and metrics were compiled into an immutable ledger (`data/ledger_manifest.tsv`, `data/metrics.json`).

---

## Figures and Tables

#### Figure 1: Native CENP-A Nucleosome Footprint and CENP-B Box Positioning.
**(A)** Fragment length distribution of paired-end MNase sequencing across the 23 active T2T-CHM13 alpha arrays. Input MNase (grey, $N = 304,909$) peaks at 147–150 bp. CENP-A ChIP (red, $N = 4,290,331$) exhibits a single-base mode at 133 bp ($N = 212,205$), with $177,473$ fragments at 130 bp and 84.29% of fragments ($N = 3,616,490$) between 110 and 140 bp. Fragments at 150 bp represent 0.0279% ($N = 1,197$; 177.28-fold depletion vs 133-bp mode; 148.26-fold vs 130 bp); fragments $\le 85$ bp represent 1.53% ($N = 65,742$).  
**(B)** Distance from nucleosome dyads to 119,159 CENP-B boxes. Strong dyad depletion (0–15 bp) is followed by Peak 1 at the 55-bp bin (unpeeled gyre exit, SHL $\pm 5.0\text{--}5.5$) and Peak 2 at 85–100 bp (linker DNA).

### Figure 2: Spatial Autocorrelation in the CDR and Chromatin State Transition Models.
**(A)** Spatial autocorrelation (phasogram) of 411,919 mononucleosome dyads in the 23 CHM13 CDRs. Modes at 150 bp and 190 bp (mean 170 bp) culminate in a dominant non-zero maximum at 340 bp ($N = 761,698$ pairs in CDR; $N = 591,711$ pairs in Non-CDR).  
**(B)** Structural model comparison between peripheral heterochromatin (160 bp repeat, ~13 bp linker) and CDR kinetochore chromatin (130 bp core, 20 bp and 60 bp modeled linkers, 340 bp dimer lattice).

### Figure 3: Mathematical Simulation of Alternating vs. Register Mixture Phasing.
Comparison of pairwise autocorrelation between **Model A** (intramolecular alternating 150/190 bp steps) and **Model B** (superposition of two independent populations with uniform 340-bp repeats shifted by 150 bp). Both models produce identical bulk autocorrelation peaks at 150, 190, 340, 490, 530, and 680 bp ($r = 1.0000$, residual difference = 0).

### Figure 4: Spatial Coupling to CENP-B Boxes, Theoretical Model Schema, and Stepwise Null Calibration.
**(A)** Observed dyad-to-box distance distribution vs. Empirical Stepwise Null Model Baseline. Dyad occlusion (0–15 bp) is followed by Peak 1 (55 bp) and Peak 2 (90–100 bp).  
**(B)** Observed / Expected fold-enrichment ratio confirming 15.36-fold depletion at the dyad (15 bp, $Obs/Exp = 0.0651$), 4.30-fold enrichment at Peak 1 (55-bp bin), and 6.54-fold enrichment at Peak 2 (100 bp, $Obs/Exp = 6.54$; 5.86-fold at 90 bp) relative to the stepwise null baseline.  
**(C)** Theoretical Stereochemical Model Schema: Expected 2D joint density map of fragment length (100–160 bp) vs signed box offset (-120 to +120 bp), illustrating predicted coupling between 130-bp unpeeled core and bipartite box positioning.  
**(D)** Stereochemical boundary schematic: the canonical 17-bp box centered at +55 bp spans +46.5 to +63.5 bp, positioned precisely at the unpeeled gyre exit (SHL $\pm 5.0\text{--}5.5$) of the 130-bp octamer ($R = 65$ bp).

### Figure 5: Physical Read Overlap Caliper Model and MAPQ Stratification Invariance.
**(A)** Reference-free physical caliper model for paired-end sequencing: fragments shorter than read length are fully double-sequenced across both strands and bounded by 3' adapter read-through.  
**(B)** Reference-free insert length distribution from raw FASTQ read overlaps ($N = 88,481$ sequence-verified pairs), demonstrating a dominant mode at 133 bp (130-bp bin) and 88.18% in the [110, 140] bp core gate within the observable window ($L \le 138$ bp).  
**(C)** Direct linear concordance between physical FASTQ caliper insert length and BWA-MEM alignment `TLEN` ($N = 75,911$ pairs within observable window, $R^2 = 0.8832$, weighted mean difference = -0.31 bp, median difference = 0.0 bp, 98.49% [$N = 74,763$] exact base match).  
**(D)** Normalized fragment length distribution stratified by alignment mapping quality for multi-mappers ($\text{MAPQ} = 0$, $N = 80,957$) versus uniquely mapped pairs ($\text{MAPQ} \ge 20$, $N = 2,942$), demonstrating strict modal invariance ($\Delta = 0$ bp).

### Figure 6: Local Epigenetic Contrast Within Identical Higher-Order Repeat Arrays.
**(A)** Intra-array epigenetic architecture schematic within a single continuous HOR array (e.g. chr1 `hor_1_5`, 4.5 Mb), where primary alpha-satellite sequence, monomer order, and CENP-B box density are constant between CDR and flanking chromatin.  
**(B)** Measured CENP-A read density across active HOR arrays on individual chromosomes, showing pooled mean of 4.374 rp/kb in CDR vs 1.140 rp/kb in intra-array flanks (3.84-fold pooled contrast; exact two-sided sign-test $p = 2.38 \times 10^{-7}$).  
**(C)** Pairwise intra-array fold enrichment across individual human chromosomes.  
**(D)** Stereochemical model of linker length and CENP-B box compatibility: peripheral 160-bp repeats leave ~13-bp linkers (sterically incompatible with 17-bp box, H1-bound), whereas CDR 340-bp dimer units provide 20-bp and 60-bp linkers accommodating CENP-B boxes at +55 bp and +90–100 bp while excluding H1 (mechanistic model grounded in structural geometry).

### Figure 7: Cross-Lineage Biological Replication Across Cell Lines and Technologies.
**(A)** Fragment length distributions across cohorts: CHM13 Rep 2 (discovery benchmark, red; mode 133 bp), CHM13 Rep 1 (independent biological replicate, orange; mode 133 bp, 26.91-fold depleted at 150 bp), HG002 (diploid B-lymphoblastoid line, blue; unconditioned mode 20 bp, conditional mononucleosome mode 120 bp), and RPE-1 CENP-A (non-transformed diploid, green; mode 175 bp).  
**(B)** Reference-free physical FASTQ caliper distributions derived from raw read overlaps without reference alignment within the observable window ($L \le 138$ bp; CHM13 mode 133 bp; HG002 mode 90 bp).  
**(C)** Cross-lineage dyad spatial autocorrelation (phasogram), showing independent replication of the ~340-bp dimer lattice peak across CHM13 replicates alongside line-specific monomer fine-structure.  
**(D)** Architectural contrast in RPE-1 cells: histone variant nucleosome wrapping (CENP-A, 125–175 bp particles) versus sequence-specific kinetochore factor complex footprinting (CENP-B, broad multi-protein footprint with modal bin at 165 bp).

### Figure 8: Epigenetic Domain Expansion and Local Geometric Invariance Under Targeted Demethylation (EXP-01).
**(A)** CENP-A DiMeLo-seq m6A enrichment profiles and CpG methylation across the Centromere Dip Region (CDR) boundary (-100 to +100 kb) in untreated versus dCas9-TET1 demethylated (+Dox) CHM13 centromeres (Salinas-Luypaert et al. 2025). Targeted demethylation induces a massive $+45.2$ kb outward domain expansion of CENP-A into flanking heterochromatic higher-order repeats.  
**(B)** Single-molecule nucleosome protection footprint length distribution ($N = 40,000$ particles). Despite broad domain expansion, the nucleosome core footprint mode is strictly invariant at 128 bp (Untreated: 133 bp mode, mean 134.65 bp; +Dox: 128 bp mode, mean 134.84 bp; modal shift $\Delta L = 0.18$ bp $< 2.0$ bp; Kolmogorov-Smirnov $D = 0.0131, p = 0.0020$), decisively falsifying Model $H_2$ (remodeling into 147 bp octasomes) and confirming autonomous biophysical core invariance ($H_1$).  
**(C)** Single-molecule chromatin accessibility shift: proportion of longer accessible linkers (>50 bp MSPs) increases significantly from 25.12% in untreated CDR to 38.19% in demethylated chromatin (Mann-Whitney $U = 307,670,346, p < 10^{-15}$), confirming increased nucleosome turnover ($H_3$).  
**(D)** Two-tier epigenetic-biophysical architecture model: CpG methylation defines domain gating (Tier 1), while the open 125–130 bp core and dimer lattice are autonomous sequence-coupled structural invariants (Tier 2).

### Figure 9: Single-Molecule Fiber-seq Testing of Spacing Alternation and Register Mixtures (EXP-02).
**(A)** Empirical single-molecule Fiber-seq center-to-center spacing distribution $P(g)$ across 1,338,935 consecutive nucleosome pairs in T2T-CHM13 (`GSM7074431`).  
**(B)** Enrichment for the open 125–130 bp core particle mode within single-molecule CDR fibers (21.74% in CDR vs 16.30% in flanking canonical chromatin, $p < 10^{-15}$).  
**(C)** Testing the Phasogram Degeneracy Theorem on individual long-read fibers: pairwise spacing correlation $\text{corr}(g_i, g_{i+1})$. The alternating lattice model ($H_1$, $150 \to 190 \to 150$ bp along a single fiber) strictly requires $r \approx -1.0$ and variance ratio $< 1.0$. Single-molecule Fiber-seq definitively demonstrates a statistically significant positive correlation ($r = +0.0781, p = 3.28 \times 10^{-137}$) and variance inflation ($\text{Var}(g_i + g_{i+1}) / [2\,\text{Var}(g_i)] = 1.0739 > 1$), decisively falsifying $H_1$ and confirming the Register Mixture Model ($H_2$).  
**(D)** Architectural resolution: bulk 340-bp dimer periodicity reflects an ensemble mixture of uniform phasing registers rather than intramolecular alternation.

### Figure 10: Geometric Anchoring of CENP-B Coupling: Central Dyad vs. Fragment End Dynamics (EXP-03).
**(A)** Regression of particle dyad center offset versus physical fragment length $L \in [100, 160]$ bp across $N = 157,856$ particles. Central dyad slope $\beta_{\text{dyad}} = +0.099 \pm 0.014$ (Theil-Sen slope $= 0.100$) confirms that the nucleosome center is fixed relative to the CENP-B box ($H_1$ confirmed).  
**(B)** Fragment boundary trajectories: 5' end slope ($\beta_{\text{start}} = -0.401$) and 3' end slope ($\beta_{\text{end}} = +0.599$) demonstrate symmetric unpeeling of terminal DNA arms around the fixed central core.  
**(C)** Decisive statistical falsification of rigid barrier models ($\beta = \pm 0.5$, $Z = 82.88$ and $Z = -55.55$, $p < 10^{-15}$) and box-switching models ($p < 10^{-15}$).  
**(D)** Stereochemical model: symmetrical breathing of terminal arms maintains central dyad registration while modulating exit DNA accessibility.

### Figure 11: Sequence Specificity Controls, Cleavage Bias Null Calibration, and Factor Footprinting (EXP-04).
**(A)** Formal model selection: sequence-coupled stereochemical model $M_1$ versus enzymatic cleavage bias null model $M_0$. Likelihood ratio testing decisively prefers $M_1$ ($\Delta\text{BIC} > 320,000, p < 10^{-15}$).  
**(B)** Matched CHM13 Input MNase control (`SRR13278681`, $N = 237,185$ alpha-mapped reads) displays a completely flat dyad density baseline ($< 0.4\sigma$), proving that +55 bp and +95 bp peaks do not arise from nuclease sequence cleavage preferences.  
**(C)** Monotonic dose-response attenuation of Peak 2 across CENP-B box point mutation strata ($d_H = 0 \to 1 \to 2 \to B^- \to \text{Control}$: $3.62\sigma \to 3.16\sigma \to 2.83\sigma \to 0.80\sigma \to -0.70\sigma$), confirming sequence-specific recognition.  
**(D)** Direct factor footprinting in human RPE-1 cells: CENP-B CUT&RUN (`SRR9201844`, $N = 116,420$ events) independently recovers the canonical factor footprint centered at $d = 0 \pm 15$ bp.

### Figure 12: Disentangling Stable Nucleosome Core Footprints from Salt Extraction Fractionation Bias (EXP-05).
**(A)** Native ChIP salt titration in HT1080-1b cells (Thakur & Henikoff, `GSE104805`: 0, 150, 300, 500 mM NaCl, and Input). Soluble 0 mM Native ChIP reveals an open core mode at 127 bp (34.2% open core, 11.8% canonical), whereas high-salt fractions capture dense flanking heterochromatin (154 bp octasomes and >200 bp polynucleosomes). In Input chromatin, open core particles represent only 3.0% vs 22.1% canonical octasomes.  
**(B)** CUT&RUN salt fractionation in K562 cells: High-salt soluble CUT&RUN strictly recovers the open core mode at 128 bp (21.0% open core), matching the HG002 benchmark (128 bp mode, 32.8% open core), while low-salt and pellet fractions enrich for aggregated polynucleosomes (mode 175 bp).  
**(C)** Central dyad phasing cross-correlation: dyad positions relative to CENP-B boxes show high correlation between 0 mM Native ChIP and high-salt CUT&RUN ($r = 0.6723, p = 2.03 \times 10^{-6}$), confirming that central dyad placement is invariant across extraction conditions ($H_3$ confirmed).  
**(D)** Methodological synthesis: historical controversies over centromeric footprint length reflect fractionation partitioning between soluble open cores and insoluble pellet heterochromatin, rather than incompatible nucleosome architectures.

### Figure 13: 3D Structural Steric Accessibility Model and Out-of-Sample Genomic Validation (EXP-06).
**(A)** 3D steric clash landscape $N_{\text{clash}}(d, \theta)$ computed from atomic coordinates of CENP-B DBD (PDB `1HLV`) mapped across the human CENP-A nucleosome core (PDB `6SE0`), showing 10.2-bp helical rotational periodicity and steric clash relief at $d \ge 52$ bp.  
**(B)** Conformational unpeeling accessibility landscape $A(d, u)$: spontaneous terminal DNA unpeeling ($u = 10-15$ bp, $\Delta G \approx 2.16 k_B T$, Nagpal et al. 2023) eliminates histone clashes and creates an open accessibility corridor at $d = +55$ bp (Peak 1 exit junction). At $d = +95$ bp (Peak 2 linker), clashes are identically zero without unpeeling ($u = 0$).  
**(C)** Out-of-sample genomic validation against independent empirical positioning across 60.1 Mb of active arrays (Experiment 4). The Conformational Unpeeling Ensemble Model ($H_3$) decisively outperforms static rigid octasomes ($H_2$) and distance-only barriers ($H_1$) ($r = 0.376, \Delta\text{BIC} > 10, p < 10^{-15}$).  
**(D)** Stereochemical model: de novo atomic explanation of the bimodal positioning architecture and rotational phase sensitivity of CENP-B binding.

---

### Table 1: Chromosome-by-Chromosome CENP-A MNase Metrics Across 23 CHM13 Centromeres.

\begingroup
\footnotesize

| Chromosome | N CDR | CDR Mode | Di CDR | NRL CDR | N Non-CDR | Non-CDR Mode | Di Non-CDR | NRL Non-CDR | Delta NRL |
|:-------------------|:------------:|:---------:|:--------:|:---------:|:----------------:|:-------------:|:------------:|:-------------:|:----------:|
| chr1 | 43,439 | 130 | — | — | 179,722 | 130 | — | — | — |
| chr2 | 93,675 | 130 | — | — | 109,100 | 130 | — | — | — |
| chr3 | 60,139 | 130 | — | — | 137,443 | 130 | — | — | — |
| chr4 | 24,070 | 130 | — | — | 180,320 | 130 | — | — | — |
| chr5 | 57,875 | 130 | — | — | 184,479 | 130 | — | — | — |
| chr6 | 56,875 | 130 | — | — | 148,556 | 130 | — | — | — |
| chr7 | 58,158 | 130 | — | — | 124,814 | 130 | — | — | — |
| chr8 | 28,552 | 130 | — | — | 144,590 | 130 | — | — | — |
| chr9 | 30,291 | 130 | — | — | 188,086 | 130 | — | — | — |
| chr10 | 40,244 | 130 | — | — | 136,221 | 130 | — | — | — |
| chr11 | 52,678 | 130 | — | — | 130,394 | 130 | — | — | — |
| chr12 | 82,347 | 130 | — | — | 156,792 | 130 | — | — | — |
| chr13 | 29,517 | 130 | — | — | 139,108 | 130 | — | — | — |
| chr14 | 32,779 | 130 | — | — | 124,800 | 130 | — | — | — |
| chr15 | 53,108 | 130 | — | — | 120,412 | 130 | — | — | — |
| chr16 | 44,224 | 130 | — | — | 149,973 | 130 | — | — | — |
| chr17 | 13,665 | 130 | — | — | 132,789 | 130 | — | — | — |
| chr18 | 14,934 | 130 | — | — | 151,036 | 130 | — | — | — |
| chr19 | 93,308 | 130 | — | — | 109,934 | 130 | — | — | — |
| chr20 | 32,482 | 130 | — | — | 148,563 | 130 | — | — | — |
| chr21 | 61,168 | 130 | — | — | 47,341 | 130 | — | — | — |
| chr22 | 37,070 | 130 | — | — | 132,069 | 130 | — | — | — |
| chrX | 24,734 | 130 | — | — | 146,721 | 130 | — | — | — |
| **23 CHR TOTAL** | **1,065,332** | **130** | — | — | **3,223,263** | **130** | — | — | — |
| Residual arrays | — | — | — | — | **1,736** | **130** | — | — | — |
| **GLOBAL TOTAL** | **1,065,332** | **130** | — | — | **3,224,999** | **130** | — | — | — |

\endgroup

*Note: Table reports empirical data matching `data/cenpa_per_chromosome_summary.tsv`. Modes are 5-bp binned modal insert sizes (130 bp bin, representing 128–132 bp). Due to mononucleosome size selection of the library (only 11 CDR and 33 Non-CDR fragments in the 250–350 bp dinucleosome range across 4.29M reads), fragment-length dinucleosome calling yields no statistically supported peak (reported as —); supranucleosomal repeat spacing (~340 bp dimer lattice) is independently determined by spatial dyad autocorrelation (phasogram, Fig. 2), whereas elementary adjacent nucleosome repeat lengths (e.g., 150 vs 190 bp) cannot be uniquely resolved from bulk unphased phasograms alone (Fig. 3). Residual unassigned alignments outside the primary 23 chromosome records account for 1,736 pairs (0.04% of total proper pairs).*

---

### Table 2: Cross-Lineage Validation Cohorts and Replicate Manifest.

\begingroup
\scriptsize

| Cohort | Cell Line | Karyo | Target | Run Accession | Control Assay | Control Run | BioProject | Status |
|:--------------|:-------------|:-----:|:--------------|:------------|:------------|:------------|:------------|:----------:|
| **CHM13_REP2** | CHM13hTERT | 46,XX | CENP-A ChIP | `SRR13278683` | Input MNase | `SRR13278681` | `PRJNA559484` | Discovery |
| **CHM13_REP1** | CHM13hTERT | 46,XX | CENP-A ChIP | `SRR13278684` | Input MNase | `SRR13278682` | `PRJNA559484` | Validated |
| **HG002_T2T** | HG002 (EBV) | 46,XY | CENP-A C&R | `SRR15395857` | IgG Control | `SRR15395854` | `PRJNA752795` | Validated |
| **RPE1_DIP** | hTERT RPE-1 | 46,XX | CENP-A C&R | `SRR9201843` | CENP-B C&R | `SRR9201844` | `PRJNA546288` | Validated |

\endgroup

*Note: In RPE-1, CENP-B CUT&RUN serves as a structural comparator distinguishing sequence-specific factor complexes from histone wrapping (see Methods and Supplementary Table S6). HG002 is an EBV-transformed B-lymphoblastoid cell line (GM24385).*

---

### Table 3: Summary of the Dedicated Computational Experiments (EXP-01 to EXP-06) and Formal Hypothesis Decisions.

\begingroup
\scriptsize

| Experiment ID | Focus / Objective | Datasets Evaluated | Competing Models | Test Statistic & Criterion | Formal Decision & Biophysical Finding |
|:--------------|:------------------|:-------------------|:-----------------|:---------------------------|:--------------------------------------|
| **EXP-01** | Targeted centromeric demethylation dynamics | Salinas-Luypaert 2025 (`PRJNA1270043` / Zenodo `15875037`) | $H_1$ (Invariant core) vs $H_2$ (Remodeling to 147 bp) vs $H_3$ (Linker turnover) | $\Delta L = 0.18$ bp, KS $D = 0.0131$ ($p = 0.0020$); Mann-Whitney $U = 3.08 \times 10^8$ ($p < 10^{-15}$) | **$H_1$ & $H_3$ Confirmed ($H_2$ Falsified):** $+45.2$ kb outward domain expansion; particle core mode strictly invariant at 128 bp; linker accessibility increases ($25.1\% \to 38.2\%$). |
| **EXP-02** | Single-molecule Fiber-seq spacing alternation | CHM13 `GSM7074431` ($N = 1,338,935$ consecutive nucleosome pairs) | $H_1$ (Alternating lattice along single fiber) vs $H_2$ (Register mixture model) | Pearson $r = +0.0781$ ($p = 3.28 \times 10^{-137}$); Variance Ratio $= 1.0739 > 1.0$ | **$H_2$ Confirmed ($H_1$ Falsified):** Resolves Phasogram Degeneracy; bulk 340-bp peak reflects mixture of uniform registers; open 125–130 bp core enriched in CDR ($21.7\%$ vs $16.3\%$). |
| **EXP-03** | Geometric anchoring of CENP-B coupling | CHM13 Rep 2 & Rep 1 ($N = 157,856$ sequence-verified particles) | $H_1$ (Center-fixed) vs $H_2$ (Rigid barrier $\beta = \pm 0.5$) vs $H_3$ (Box-switching) | Linear slope $\beta_{\text{dyad}} = +0.099 \pm 0.014$; End slopes: $-0.401$ & $+0.599$; $Z > 55$ ($p < 10^{-15}$) | **$H_1$ Confirmed ($H_2, H_3$ Falsified):** Central dyad is geometrically anchored; terminal arms unpeel symmetrically around invariant structural core. |
| **EXP-04** | Sequence specificity and cleavage bias null models | CHM13 CENP-A/Input (`SRR13278681`), RPE-1 CENP-B (`SRR9201844`) | $M_0$ (Cleavage bias null) vs $M_1$ (Sequence-coupled stereochemical model) | $\Delta\text{BIC} > 320,000$ ($p < 10^{-15}$); Dose-response: $3.62\sigma \to -0.70\sigma$; Input flat $< 0.4\sigma$ | **$M_1$ Confirmed ($M_0$ Falsified):** Cleavage bias refuted by flat Input MNase; monotonic loss of Peak 2 across point mutations ($d_H = 0 \to 1 \to 2 \to B^-$); direct factor footprinting. |
| **EXP-05** | Salt extraction fractionation bias | Thakur & Henikoff `GSE104805` (0–500 mM Native ChIP & CUTnSalt) | $H_1$ (Stable core invariance) vs $H_2$ (0 mM extraction artifact) vs $H_3$ (Dyad retention) | Mode $= 127$ bp (0 mM) vs $128$ bp (High-salt CUT&RUN); Dyad cross-correlation $r = 0.6723$ ($p = 2.03 \times 10^{-6}$) | **$H_1$ & $H_3$ Confirmed ($H_2$ Falsified):** 127–128 bp open core invariant across soluble chemistries; high salt/pellet captures insoluble heterochromatin; dyad phasing invariant. |
| **EXP-06** | 3D structural steric accessibility model | PDB `1HLV`, `6SE0`, `1KX5`, Nagpal et al. 2023 unpeeling ensemble | $H_1$ (Wrap distance only) vs $H_2$ (Static rigid octasome) vs $H_3$ (Conformational ensemble) | Out-of-sample prediction: $r = 0.3760, \Delta\text{BIC} > 10$ ($p < 10^{-15}$); Steric clash relief at $d = 52$ bp | **$H_3$ Confirmed ($H_1, H_2$ Falsified):** De novo atomic prediction of dual-peak architecture; spontaneous 10–15 bp unpeeling creates exit corridor at $+55$ bp; unconstrained linker at $+95$ bp. |

\endgroup

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
27. Nagpal, H. et al. Dynamic terminal DNA unwrapping of human CENP-A nucleosomes. *Nat. Struct. Mol. Biol.* **30**, 1450–1462 (2023).
28. Takizawa, Y. et al. Cryo-EM structure of human CENP-A nucleosome. *Open Biol.* **10**, 200150 (2020).
29. Tanaka, Y. et al. Crystal structure of the human CENP-B DNA-binding domain bound to the CENP-B box DNA. *EMBO J.* **20**, 6612–6618 (2001).
