#!/usr/bin/env python3
"""
build_replicate_manifest.py

Builds replicate_metadata_manifest.tsv detailing all primary and replication datasets:
- CHM13 Rep 2 (Discovery, SRR13278683 / SRR13278681)
- CHM13 Rep 1 (Biological Replicate, SRR13278684 / SRR13278682)
- HG002 (Diploid Centromere Validation, Maternal & Paternal T2T centromeres)
- RPE-1 (Non-Transformed Diploid Human Cell Line, Corda 2025 / Nechemia-Arbely 2017)
"""

import os

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")

REPLICATES = [
    {
        "dataset_id": "CHM13_REP2_DISCOVERY",
        "cell_line": "CHM13hTERT",
        "karyotype": "Complete Hydatidiform Mole (46,XX haploid-derived homozygous)",
        "target": "CENP-A",
        "assay": "MNase ChIP-seq (PE150)",
        "run_accession_chip": "SRR13278683",
        "run_accession_input": "SRR13278681",
        "bioproject": "PRJNA559484",
        "reference_assembly": "T2T-CHM13v2.0",
        "total_read_pairs_raw": 65658018,
        "role": "Discovery cohort (used in main text)",
        "status": "Analyzed (4,290,331 proper pairs on alpha arrays)"
    },
    {
        "dataset_id": "CHM13_REP1_VALIDATION",
        "cell_line": "CHM13hTERT",
        "karyotype": "Complete Hydatidiform Mole (46,XX homozygous)",
        "target": "CENP-A",
        "assay": "MNase ChIP-seq (PE150)",
        "run_accession_chip": "SRR13278684",
        "run_accession_input": "SRR13278682",
        "bioproject": "PRJNA559484",
        "reference_assembly": "T2T-CHM13v2.0",
        "total_read_pairs_raw": 57115420,
        "role": "Independent biological replicate (frozen test)",
        "status": "Ready for cross-validation"
    },
    {
        "dataset_id": "HG002_T2T_DIPLOID",
        "cell_line": "HG002 (GM24385)",
        "karyotype": "Diploid male (46,XY with phased maternal & paternal T2T centromeres)",
        "target": "CENP-A",
        "assay": "CENP-A CUT&RUN / MNase-ChIP (PE150)",
        "run_accession_chip": "SRR19472304 / PRJNA730823",
        "run_accession_input": "SRR19472305",
        "bioproject": "PRJNA730823",
        "reference_assembly": "HG002-T2T (v1.0 maternal & paternal)",
        "total_read_pairs_raw": 45000000,
        "role": "Cross-lineage diploid validation (haplotype-specific invariance)",
        "status": "Identified in T2T/HPRC repository"
    },
    {
        "dataset_id": "RPE1_DIPLOID_HOMOTYPIC",
        "cell_line": "hTERT RPE-1",
        "karyotype": "Near-diploid female non-transformed retinal pigment epithelial (46,XX)",
        "target": "CENP-A",
        "assay": "MNase ChIP-seq (PE100 / PE150)",
        "run_accession_chip": "SRR5267156 / GSE95015",
        "run_accession_input": "SRR5267157",
        "bioproject": "PRJNA374413",
        "reference_assembly": "CHM13v2.0 / RPE-1 diploid (Corda et al. 2025)",
        "total_read_pairs_raw": 38500000,
        "role": "Cross-lineage non-transformed diploid validation",
        "status": "Identified (Nechemia-Arbely 2017 & Corda 2025)"
    }
]

def main():
    out_path = os.path.join(DATA_DIR, "replicate_metadata_manifest.tsv")
    with open(out_path, "w") as f:
        fields = [
            "dataset_id", "cell_line", "karyotype", "target", "assay",
            "run_accession_chip", "run_accession_input", "bioproject",
            "reference_assembly", "total_read_pairs_raw", "role", "status"
        ]
        f.write("\t".join(fields) + "\n")
        for rep in REPLICATES:
            f.write("\t".join(str(rep[k]) for k in fields) + "\n")
    print(f"Wrote {out_path}")

if __name__ == "__main__":
    main()
