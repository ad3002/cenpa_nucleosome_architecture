#!/usr/bin/env bash
# 01_fetch_and_align.sh
# Downloads paired-end MNase sequencing slices for CHM13 Input (SRR13278681)
# and CENP-A ChIP-seq (SRR13278683), aligns via BWA-MEM, and indexes.

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"
DATA_DIR="$BASE_DIR/data"
RAW_DIR="$BASE_DIR/raw_cache"

mkdir -p "$RAW_DIR" "$DATA_DIR"

THREADS=${THREADS:-16}
PAIRS_LIMIT=${PAIRS_LIMIT:-3000000}

# 1. References
ALPHA_FA="$RAW_DIR/chm13_alpha_arrays.fa"
if [ ! -f "$ALPHA_FA" ]; then
    echo "Downloading / extracting CHM13 alpha-satellite arrays reference..."
    # Can be extracted from CHM13v2.0 2bit or downloaded
    if [ -f "/mnt/data/claude/2026-09-14_nucleosome_genomics/cache/chm13_alpha_arrays.fa" ]; then
        cp "/mnt/data/claude/2026-09-14_nucleosome_genomics/cache/chm13_alpha_arrays.fa" "$ALPHA_FA"
    fi
fi

if [ -f "$ALPHA_FA" ] && [ ! -f "$ALPHA_FA.bwt" ]; then
    echo "Building BWA index on $ALPHA_FA..."
    bwa index "$ALPHA_FA"
fi

align_sra() {
    local ACC=$1
    local VOL=$2
    local OUT_BAM="$RAW_DIR/${ACC}.sorted.bam"

    if [ -s "$OUT_BAM" ]; then
        echo "BAM for $ACC already exists: $OUT_BAM"
        return
    fi

    local U1="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR132/${VOL}/${ACC}/${ACC}_1.fastq.gz"
    local U2="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR132/${VOL}/${ACC}/${ACC}_2.fastq.gz"

    local R1="$RAW_DIR/${ACC}_raw_1.fastq.gz"
    local R2="$RAW_DIR/${ACC}_raw_2.fastq.gz"
    local S1="$RAW_DIR/${ACC}_sync_1.fastq.gz"
    local S2="$RAW_DIR/${ACC}_sync_2.fastq.gz"

    echo "Downloading slice for $ACC (600 MB)..."
    [ ! -s "$R1" ] && curl -4 -s -r 0-629145600 --retry 5 "$U1" -o "$R1"
    [ ! -s "$R2" ] && curl -4 -s -r 0-629145600 --retry 5 "$U2" -o "$R2"

    echo "Synchronizing paired-end records for $ACC..."
    python3 "$SCRIPT_DIR/sync_paired.py" "$R1" "$R2" "$S1" "$S2" "$PAIRS_LIMIT"

    echo "Aligning $ACC with BWA-MEM ($THREADS threads)..."
    bwa mem -t "$THREADS" "$ALPHA_FA" "$S1" "$S2" \
        | samtools view -b -F 4 - \
        | samtools sort -@ "$THREADS" -o "$OUT_BAM" -
    samtools index "$OUT_BAM"
    echo "Completed $ACC alignment."
}

echo "=== Processing Input MNase (SRR13278681) ==="
align_sra "SRR13278681" "081"

echo "=== Processing CENP-A ChIP-seq (SRR13278683) ==="
align_sra "SRR13278683" "083"

echo "Raw data alignment complete."
