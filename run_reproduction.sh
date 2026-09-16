#!/usr/bin/env bash
# run_reproduction.sh
# Master reproduction script for "Human CENP-A Nucleosomes Form an Open 125-130 bp Particle
# Phased to a 340 bp Alpha-Satellite Dimer Lattice with Bipartite CENP-B Linker Geometry".

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODE="quick"

if [ "$1" == "--full-raw" ]; then
    MODE="full"
fi

echo "================================================================================"
echo "  CENP-A Nucleosome Architecture — Reproduction Suite"
echo "  Mode: $MODE"
echo "================================================================================"

if [ "$MODE" == "full" ]; then
    echo "1. Fetching raw SRA slices and running BWA-MEM alignments..."
    bash "$SCRIPT_DIR/scripts/01_fetch_and_align.sh"

    echo "2. Analyzing particle geometries from BAM files..."
    python3 "$SCRIPT_DIR/scripts/02_analyze_particles.py" \
        "$SCRIPT_DIR/raw_cache/SRR13278681.sorted.bam" \
        "$SCRIPT_DIR/data/chm13_cdr_intervals.bed" \
        "$SCRIPT_DIR/data/chm13_cenpb_boxes_coords.tsv" \
        "$SCRIPT_DIR/data/input"

    python3 "$SCRIPT_DIR/scripts/02_analyze_particles.py" \
        "$SCRIPT_DIR/raw_cache/SRR13278683.sorted.bam" \
        "$SCRIPT_DIR/data/chm13_cdr_intervals.bed" \
        "$SCRIPT_DIR/data/chm13_cenpb_boxes_coords.tsv" \
        "$SCRIPT_DIR/data/cenpa"

    echo "3. Computing spatial phasograms..."
    python3 "$SCRIPT_DIR/scripts/03_compute_phasogram.py" \
        "$SCRIPT_DIR/raw_cache/SRR13278683.sorted.bam" \
        "$SCRIPT_DIR/data/chm13_cdr_intervals.bed" \
        "$SCRIPT_DIR/data/cenpa_cdr_phasogram.tsv"
else
    echo "1. Quick mode: Using curated data tables in data/ (instant reproduction)..."
    echo "   (To run full raw alignment from EBI/SRA, pass: ./run_reproduction.sh --full-raw)"
fi

echo "Generating publication figures (PNG, PDF, SVG)..."
python3 "$SCRIPT_DIR/scripts/04_plot_figures.py"

echo "================================================================================"
echo "  Reproduction complete!"
echo "  Figures saved to: paper/figures/"
echo "    - Fig1_cenpa_core_and_box_geometry.png (.pdf, .svg)"
echo "    - Fig2_cdr_phasogram_and_chromatin_state.png (.pdf, .svg)"
echo "  Manuscript available at: paper/manuscript.md"
echo "================================================================================"
