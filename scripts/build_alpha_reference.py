#!/usr/bin/env python3
"""
build_alpha_reference.py

Extracts active centromeric HOR arrays for all 23 human chromosomes from T2T-CHM13v2.0 (hs1)
using UCSC twoBitToFa, and produces raw_cache/chm13_alpha_arrays.fa formatted with
the accession-coordinate headers required by 02_analyze_particles.py.
"""

import os
import sys
import subprocess
import urllib.request

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(SCRIPT_DIR)
RAW_DIR = os.path.join(BASE_DIR, "raw_cache")
DATA_DIR = os.path.join(BASE_DIR, "data")
os.makedirs(RAW_DIR, exist_ok=True)

TWOBIT_BIN = os.path.join(RAW_DIR, "twoBitToFa")
HS1_2BIT_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.2bit"
CENSAT_URL = "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.0.bed"
CDR_BED = os.path.join(DATA_DIR, "chm13_cdr_intervals.bed")
OUT_FA = os.path.join(RAW_DIR, "chm13_alpha_arrays.fa")

# Chromosome mapping
chrom_map = {f"chr{i}": f"NC_{60925 + i - 1:06d}.1" for i in range(1, 23)}
chrom_map["chrX"] = "NC_060947.1"

def main():
    if os.path.exists(OUT_FA) and os.path.getsize(OUT_FA) > 10000000:
        print(f"Reference {OUT_FA} already exists ({os.path.getsize(OUT_FA):,} bytes).")
        return

    # 1. Load CDR intervals
    cdrs = {}
    with open(CDR_BED) as f:
        for line in f:
            p = line.strip().split()
            if len(p) >= 3:
                chrom = p[0].replace("_hap1", "")
                cdrs[chrom] = (int(p[1]), int(p[2]))

    # 2. Download CenSat BED and identify active HOR per chromosome
    print("Fetching CenSat annotations to identify active HOR arrays...")
    active_hors = {}
    req = urllib.request.urlopen(CENSAT_URL)
    for raw_line in req:
        line = raw_line.decode("utf-8")
        if line.startswith("track") or not line.strip():
            continue
        p = line.strip().split("\t")
        if len(p) < 4:
            continue
        c, s, e, name = p[0], int(p[1]), int(p[2]), p[3]
        if c in cdrs:
            cs, ce = cdrs[c]
            if max(s, cs) < min(e, ce):
                active_hors[c] = (s, e, name)

    print(f"Found {len(active_hors)} active HOR arrays overlapping CDRs.")

    # 3. Extract sequences using twoBitToFa
    temp_fa = OUT_FA + ".tmp"
    with open(temp_fa, "w") as out_f:
        for c in sorted(active_hors.keys(), key=lambda x: (int(x[3:]) if x[3:].isdigit() else 99)):
            s, e, name = active_hors[c]
            acc = chrom_map[c]
            header = f"{acc}_{s}_{e}"
            print(f"Extracting {c} ({acc}) {s:,} - {e:,} ({name})...")
            
            cmd = [TWOBIT_BIN, HS1_2BIT_URL, f"-seq={c}", f"-start={s}", f"-end={e}", "stdout"]
            proc = subprocess.run(cmd, stdout=subprocess.PIPE, text=True, check=True)
            
            lines = proc.stdout.strip().split("\n")
            out_f.write(f">{header}\n")
            for seq_line in lines[1:]:
                out_f.write(seq_line + "\n")

    os.replace(temp_fa, OUT_FA)
    print(f"Successfully generated {OUT_FA} ({os.path.getsize(OUT_FA):,} bytes).")

    # 4. Build BWA index
    print("Building BWA index...")
    subprocess.run(["bwa", "index", OUT_FA], check=True)
    print("BWA index build complete.")

if __name__ == "__main__":
    main()
