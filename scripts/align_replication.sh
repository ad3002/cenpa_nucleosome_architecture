#!/usr/bin/env bash
set -eo pipefail
ALPHA_FA="raw_cache/chm13_alpha_arrays.fa"
REP_DIR="raw_cache/replication"
THREADS=8

for ACC in SRR13278684 SRR15395857 SRR9201843 SRR9201844; do
    S1="$REP_DIR/${ACC}_sync_1.fastq.gz"
    S2="$REP_DIR/${ACC}_sync_2.fastq.gz"
    OUT_BAM="$REP_DIR/${ACC}_slice.sorted.bam"
    
    if [ -s "$OUT_BAM" ]; then
        echo "$OUT_BAM already exists."
        continue
    fi
    
    echo "Aligning $ACC to $ALPHA_FA with $THREADS threads..."
    bwa mem -t "$THREADS" "$ALPHA_FA" "$S1" "$S2" 2>/dev/null \
        | samtools view -b -F 4 - \
        | samtools sort -@ "$THREADS" -o "$OUT_BAM" -
    samtools index "$OUT_BAM"
    echo "Completed $ACC: $(samtools view -c "$OUT_BAM") mapped records."
done
echo "All replication alignments complete."
