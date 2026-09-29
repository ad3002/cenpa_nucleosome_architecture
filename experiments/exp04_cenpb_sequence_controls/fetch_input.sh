#!/usr/bin/env bash
# fetch_input.sh - Fetch matched CHM13 Input MNase (SRR13278681) control and align to alpha arrays
set -eo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
RAW_DIR="$BASE_DIR/raw_cache"
SCRIPT_DIR="$BASE_DIR/scripts"
ALPHA_FA="$RAW_DIR/chm13_alpha_arrays.fa"
THREADS=16
ACC="SRR13278681"
VOL="081"

OUT_BAM="$RAW_DIR/${ACC}_slice.sorted.bam"

if [ -s "$OUT_BAM" ]; then
    echo "$OUT_BAM already exists. Mapped count: $(samtools view -c "$OUT_BAM")"
    exit 0
fi

echo "Fetching slice of $ACC from ENA..."
R1="$RAW_DIR/${ACC}_raw_1.fastq.gz"
R2="$RAW_DIR/${ACC}_raw_2.fastq.gz"
S1="$RAW_DIR/${ACC}_sync_1.fastq.gz"
S2="$RAW_DIR/${ACC}_sync_2.fastq.gz"

U1="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR132/${VOL}/${ACC}/${ACC}_1.fastq.gz"
U2="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR132/${VOL}/${ACC}/${ACC}_2.fastq.gz"

# 60MB slice corresponds to ~1.5M read pairs
[ ! -s "$R1" ] && curl -4 -s -r 0-62914560 --retry 5 "$U1" -o "$R1"
[ ! -s "$R2" ] && curl -4 -s -r 0-62914560 --retry 5 "$U2" -o "$R2"

echo "Synchronizing paired records..."
python3 "$SCRIPT_DIR/sync_paired.py" "$R1" "$R2" "$S1" "$S2" 3000000

echo "Aligning $ACC with BWA-MEM ($THREADS threads)..."
bwa mem -t "$THREADS" "$ALPHA_FA" "$S1" "$S2" \
    | samtools view -b -F 4 - \
    | samtools sort -@ "$THREADS" -o "$OUT_BAM" -
samtools index "$OUT_BAM"

echo "Alignment complete for $ACC: $(samtools view -c "$OUT_BAM") records mapped."
