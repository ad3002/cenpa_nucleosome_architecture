#!/usr/bin/env bash
# run_experiment.sh - Experiment 4 runner
set -eo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"
REPO_ROOT="$(dirname "$BASE_DIR")"

echo "=== Running Experiment 4: CENP-B Sequence Controls & Cleavage Bias ==="
if [ -d "$REPO_ROOT/.venv" ]; then
    "$REPO_ROOT/.venv/bin/python" "$SCRIPT_DIR/run_exp04.py"
elif command -v python3 >/dev/null 2>&1; then
    python3 "$SCRIPT_DIR/run_exp04.py"
else
    echo "Error: Python 3 not found." >&2
    exit 1
fi
echo "=== Experiment 4 Complete ==="
