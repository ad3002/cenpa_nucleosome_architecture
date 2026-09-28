#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

echo "================================================================================"
echo "  Executing Experiment 3: CENP-B Coupling Geometric Anchor Analysis"
echo "  (Dyad Center vs. Fragment Termini Slopes & Monomer Box Switching)"
echo "================================================================================"

python3 "$SCRIPT_DIR/analyze_dyad_vs_ends.py"

echo "Generated artifacts in:"
echo "  - Data:    $SCRIPT_DIR/data/"
echo "  - Figures: $SCRIPT_DIR/figures/"

echo "Experiment 3 successfully executed!"
