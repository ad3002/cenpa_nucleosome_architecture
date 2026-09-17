#!/usr/bin/env python3
"""
annotate_cenpb_boxes.py

Annotates canonical 17-bp CENP-B boxes ([CT]TTCGTTGGAA[AG]CGGGA) across centromeric
alpha-satellite arrays on both forward and reverse strands.
Produces coordinate TSV for downstream dyad-to-box distance analysis.
"""

import sys
import os
import re
import argparse

# Canonical 17-bp CENP-B box motif:
# Forward: 5'-[CT]TTCGTTGGAA[AG]CGGGA-3' (YTTCGTTGGAARCGGGA)
# Reverse: 5'-TCCCG[CT]TTCCAACGAA[AG]-3' (TCCCGYTTCCAACGAAR)
FWD_REGEX = re.compile(r"[CT]TTCGTTGGAA[AG]CGGGA", re.IGNORECASE)
REV_REGEX = re.compile(r"TCCCG[CT]TTCCAACGAA[AG]", re.IGNORECASE)

def annotate_fasta(fasta_path, out_tsv):
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"Input FASTA file not found: {fasta_path}")

    os.makedirs(os.path.dirname(os.path.abspath(out_tsv)), exist_ok=True)

    total_boxes = 0
    total_arrays = 0

    with open(fasta_path, "r") as f_in, open(out_tsv, "w") as f_out:
        f_out.write("array_id\tbox_start\tbox_end\tstrand\tsequence\n")
        
        current_id = None
        current_seq = []

        def process_record(rec_id, seq_str):
            nonlocal total_boxes
            count = 0
            # Search forward strand
            for match in FWD_REGEX.finditer(seq_str):
                f_out.write(f"{rec_id}\t{match.start()}\t{match.end()}\t+\t{match.group(0).upper()}\n")
                count += 1
            # Search reverse strand
            for match in REV_REGEX.finditer(seq_str):
                f_out.write(f"{rec_id}\t{match.start()}\t{match.end()}\t-\t{match.group(0).upper()}\n")
                count += 1
            total_boxes += count

        for line in f_in:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id:
                    process_record(current_id, "".join(current_seq))
                    total_arrays += 1
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        if current_id:
            process_record(current_id, "".join(current_seq))
            total_arrays += 1

    print(f"Annotated {total_boxes:,} CENP-B boxes across {total_arrays:,} alpha-satellite arrays.")
    print(f"Saved coordinates to: {out_tsv}")

def main():
    parser = argparse.ArgumentParser(description="Annotate 17-bp CENP-B boxes in centromeric alpha arrays FASTA.")
    parser.add_argument("fasta", nargs="?", default=None, help="Input FASTA of alpha-satellite arrays")
    parser.add_argument("out_tsv", nargs="?", default=None, help="Output TSV of annotated CENP-B box coordinates")
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.dirname(script_dir)
    default_fasta = os.path.join(base_dir, "raw_cache", "chm13_alpha_arrays.fa")
    default_tsv = os.path.join(base_dir, "raw_cache", "chm13_cenpb_boxes_coords.tsv")

    fasta_path = args.fasta if args.fasta else default_fasta
    out_tsv = args.out_tsv if args.out_tsv else default_tsv

    if not os.path.exists(fasta_path):
        print(f"Notice: Target FASTA {fasta_path} not found. Ensure raw reference is fetched first.")
        sys.exit(0)

    annotate_fasta(fasta_path, out_tsv)

if __name__ == "__main__":
    main()
