#!/usr/bin/env bash
# run_experiment.sh - Experiment 2 runner
set -eo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"
REPO_ROOT="$(dirname "$BASE_DIR")"

echo "=== Running Experiment 2: Single-Molecule Fiber-seq Spacing Alternation ==="
if [ -d "$REPO_ROOT/.venv" ]; then
    "$REPO_ROOT/.venv/bin/python" "$SCRIPT_DIR/analyze_single_molecule_fiberseq.py"
elif command -v python3 >/dev/null 2>&1; then
    python3 "$SCRIPT_DIR/analyze_single_molecule_fiberseq.py"
else
    echo "Error: Python 3 not found." >&2
    exit 1
fi
echo "=== Experiment 2 Complete ==="
