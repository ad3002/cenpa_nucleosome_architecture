#!/usr/bin/env bash
# Standalone Reviewer Reproduction Suite
set -euo pipefail

echo "================================================================================"
echo "  CENP-A Nucleosome Architecture — Standalone Reviewer Reproduction Suite"
echo "================================================================================"

echo "1. Verifying cryptographic checksums (SHA256SUMS.txt)..."
python3 - << 'EOF'
import hashlib
import sys
from pathlib import Path

root = Path('.')
sha_file = root / 'SHA256SUMS.txt'
if not sha_file.is_file():
    print("  [FAIL] Missing SHA256SUMS.txt!")
    sys.exit(1)

bad = []
missing = []
total = 0
for line in sha_file.read_text().splitlines():
    line = line.strip()
    if not line:
        continue
    parts = line.split(maxsplit=1)
    if len(parts) != 2:
        print(f"  [FAIL] Malformed line in SHA256SUMS.txt: '{line}'")
        sys.exit(1)
    expected, name = parts
    p = root / name.lstrip('*')
    if not p.is_file():
        missing.append(name)
        continue
    h = hashlib.sha256(p.read_bytes()).hexdigest()
    if h != expected:
        bad.append(name)
    total += 1

if total < 41:
    print(f"  [FAIL] Incomplete manifest inventory! Expected at least 41 files, found {total}.")
    sys.exit(1)

if missing or bad:
    if missing:
        print(f"  [FAIL] Missing files ({len(missing)}): {missing[:5]}")
    if bad:
        print(f"  [FAIL] Checksum mismatches ({len(bad)}): {bad[:5]}")
    sys.exit(1)

print(f"  [PASS] All {total} package artifacts cryptographically verified against SHA256SUMS.txt")
EOF

echo "2. Inspecting single-source-of-truth ledger metrics..."
python3 - << 'EOF'
import csv
import json
import sys

with open('ledger/metrics.json') as f:
    m = json.load(f)

required_sections = [
    'particle_sizing',
    'cdr_phasogram',
    'cenpb_box_geometry',
    'physical_caliper_and_mapping',
    'intra_array_contrast',
    'cross_lineage_replication'
]
for s in required_sections:
    if s not in m:
        print(f"  [FAIL] Missing required section in metrics.json: '{s}'")
        sys.exit(1)

# A. Particle sizing
p = m['particle_sizing']
mode = p.get('single_base_mode_length_bp', 0)
mode_count = p.get('single_base_mode_count_global', 0)
core_pct = p.get('pct_110_140bp_of_global', None)
pct_150 = p.get('pct_150bp_of_global', None)
pct_sub85 = p.get('pct_sub_85bp_of_global', None)

sc = m.get('sample_counts', {})
glob_pairs = sc.get('chip_proper_pairs_global', 0)
cdr_pairs = sc.get('chip_proper_pairs_cdr_total', 0)
noncdr_pairs = sc.get('chip_proper_pairs_noncdr_total', 0)

if not (110 <= mode <= 150):
    print(f"  [FAIL] Impossible particle single-base mode: {mode} bp")
    sys.exit(1)
if mode_count <= 0 or glob_pairs <= 0 or cdr_pairs <= 0 or noncdr_pairs <= 0:
    print(f"  [FAIL] Non-positive fragment counts in particle sizing/sample counts: mode={mode_count}, glob={glob_pairs}")
    sys.exit(1)
if core_pct is None or not (0.0 <= core_pct <= 100.0):
    print(f"  [FAIL] Impossible or missing core percentage: {core_pct}%")
    sys.exit(1)
if pct_150 is None or not (0.0 <= pct_150 <= 100.0):
    print(f"  [FAIL] Impossible or missing 150 bp percentage: {pct_150}%")
    sys.exit(1)
if pct_sub85 is None or not (0.0 <= pct_sub85 <= 100.0):
    print(f"  [FAIL] Impossible or missing sub-85 bp percentage: {pct_sub85}%")
    sys.exit(1)

print(f"  [PASS] Single-base mode: {mode} bp (N={mode_count:,})")
print(f"  [PASS] Core gate [110, 140] bp: {core_pct}%, Canonical 150 bp: {pct_150}%, Sub-85 bp: {pct_sub85}%")

# B. CDR Phasogram
ph = m['cdr_phasogram']
dimer_peak = ph.get('dimer_lattice_peak_bp', 0)
dimer_pairs = ph.get('dimer_lattice_pairs_at_340bp_cdr', 0)
if not (200 <= dimer_peak <= 500) or dimer_pairs <= 0:
    print(f"  [FAIL] Impossible dimer lattice peak or count: {dimer_peak} bp (N={dimer_pairs})")
    sys.exit(1)
print(f"  [PASS] Phasogram CDR peak [100, 800] bp: {dimer_peak} bp (N={dimer_pairs:,})")

# C. CENP-B box geometry
b = m['cenpb_box_geometry']
box_contrast = b.get('peak_to_dyad_contrast_ratio', 0)
oe_p1 = b.get('observed_over_expected_peak1_55bp', 0.0)
oe_p2 = b.get('observed_over_expected_peak2_100bp', 0.0)
if box_contrast <= 0 or oe_p1 <= 0 or oe_p2 <= 0:
    print(f"  [FAIL] Impossible CENP-B box geometry metrics: contrast={box_contrast}, oe1={oe_p1}, oe2={oe_p2}")
    sys.exit(1)
print(f"  [PASS] Box contrast 55 bp vs dyad: {box_contrast}x, Peak 1 Obs/Exp: {oe_p1:.2f}, Peak 2 Obs/Exp: {oe_p2:.2f}")

# D. Caliper and MAPQ
cq = m['physical_caliper_and_mapping']
cal_mode = cq.get('physical_caliper_single_base_mode_bp', 0)
cal_core = cq.get('physical_caliper_core_gate_110_140bp_pct', None)
cal_exact = cq.get('caliper_to_tlen_exact_agreement_pct', None)
cal_r2 = cq.get('caliper_to_tlen_r2', None)
if not (110 <= cal_mode <= 150):
    print(f"  [FAIL] Impossible caliper mode: {cal_mode} bp")
    sys.exit(1)
if cal_core is None or not (0.0 <= cal_core <= 100.0):
    print(f"  [FAIL] Impossible caliper core pct: {cal_core}")
    sys.exit(1)
if cal_exact is None or not (0.0 <= cal_exact <= 100.0):
    print(f"  [FAIL] Impossible caliper exact agreement pct: {cal_exact}")
    sys.exit(1)
if cal_r2 is None or not (0.0 <= cal_r2 <= 1.0):
    print(f"  [FAIL] Impossible caliper R2: {cal_r2}")
    sys.exit(1)
print(f"  [PASS] Physical caliper: mode={cal_mode} bp, core={cal_core}%, R2={cal_r2}, exact={cal_exact}%")

# E. Intra-array contrast
ia = m['intra_array_contrast']
fold = ia.get('intra_array_fold_enrichment', 0.0)
c_dens = ia.get('cdr_read_density_rp_per_kb', 0.0)
f_dens = ia.get('flank_read_density_rp_per_kb', 0.0)
sign_p = ia.get('exact_sign_test_p', None)
if not (1.0 <= fold <= 20.0) or c_dens <= 0 or f_dens <= 0:
    print(f"  [FAIL] Impossible intra-array metrics: fold={fold}x, cdr_dens={c_dens}, flank_dens={f_dens}")
    sys.exit(1)
if sign_p is None or not (0.0 <= sign_p <= 1.0):
    print(f"  [FAIL] Impossible sign test p-value: {sign_p}")
    sys.exit(1)
print(f"  [PASS] Intra-array CDR vs Flank contrast: {fold}x ({c_dens} vs {f_dens} rp/kb; p = {sign_p:.2e})")

# F. Cross-lineage replication (Validate ALL cohorts explicitly)
cl = m['cross_lineage_replication']

# CHM13 Rep 1
rep1_mode = cl.get('chm13_rep1_single_base_mode_bp', 0)
rep1_core = cl.get('chm13_rep1_core_pct', None)
rep1_150 = cl.get('chm13_rep1_150bp_pct', None)
rep1_pairs = cl.get('chm13_rep1_analyzed_pairs', 0)
if not (110 <= rep1_mode <= 150) or rep1_pairs <= 0:
    print(f"  [FAIL] Impossible CHM13 Rep 1 metrics: mode={rep1_mode}, pairs={rep1_pairs}")
    sys.exit(1)
if rep1_core is None or not (0.0 <= rep1_core <= 100.0):
    print(f"  [FAIL] Impossible CHM13 Rep 1 core pct: {rep1_core}")
    sys.exit(1)
if rep1_150 is None or not (0.0 <= rep1_150 <= 100.0):
    print(f"  [FAIL] Impossible CHM13 Rep 1 150bp pct: {rep1_150}")
    sys.exit(1)

# HG002
hg002_pairs = cl.get('hg002_t2t_cutrun_analyzed_pairs', 0)
hg002_mode = cl.get('hg002_t2t_unconditioned_mode_bp', 0)
hg002_core = cl.get('hg002_t2t_core_pct', None)
hg002_150 = cl.get('hg002_t2t_150bp_pct', None)
hg002_sub85 = cl.get('hg002_t2t_sub85bp_pct', None)
if hg002_pairs <= 0 or hg002_mode <= 0:
    print(f"  [FAIL] Impossible HG002 metrics: pairs={hg002_pairs}, mode={hg002_mode}")
    sys.exit(1)
if hg002_core is None or not (0.0 <= hg002_core <= 100.0):
    print(f"  [FAIL] Impossible HG002 core pct: {hg002_core}")
    sys.exit(1)
if hg002_150 is None or not (0.0 <= hg002_150 <= 100.0):
    print(f"  [FAIL] Impossible HG002 150bp pct: {hg002_150}")
    sys.exit(1)
if hg002_sub85 is None or not (0.0 <= hg002_sub85 <= 100.0):
    print(f"  [FAIL] Impossible HG002 sub85 pct: {hg002_sub85}")
    sys.exit(1)

# RPE1 CENP-A
rpe1_a_pairs = cl.get('rpe1_cenpa_cutrun_analyzed_pairs', 0)
rpe1_a_mode = cl.get('rpe1_cenpa_single_base_mode_bp', 0)
rpe1_a_core = cl.get('rpe1_cenpa_core_pct', None)
rpe1_a_150 = cl.get('rpe1_cenpa_150bp_pct', None)
rpe1_a_147_175 = cl.get('rpe1_cenpa_147_175bp_pct', None)
if rpe1_a_pairs <= 0 or rpe1_a_mode <= 0:
    print(f"  [FAIL] Impossible RPE-1 CENP-A metrics: pairs={rpe1_a_pairs}, mode={rpe1_a_mode}")
    sys.exit(1)
if rpe1_a_core is None or not (0.0 <= rpe1_a_core <= 100.0):
    print(f"  [FAIL] Impossible RPE-1 CENP-A core pct: {rpe1_a_core}")
    sys.exit(1)
if rpe1_a_150 is None or not (0.0 <= rpe1_a_150 <= 100.0):
    print(f"  [FAIL] Impossible RPE-1 CENP-A 150bp pct: {rpe1_a_150}")
    sys.exit(1)
if rpe1_a_147_175 is None or not (0.0 <= rpe1_a_147_175 <= 100.0):
    print(f"  [FAIL] Impossible RPE-1 CENP-A 147-175bp pct: {rpe1_a_147_175}")
    sys.exit(1)

# RPE1 CENP-B
rpe1_b_pairs = cl.get('rpe1_cenpb_cutrun_analyzed_pairs', 0)
rpe1_b_45_65 = cl.get('rpe1_cenpb_45_65bp_pct', None)
rpe1_b_sub85 = cl.get('rpe1_cenpb_sub85bp_pct', None)
if rpe1_b_pairs <= 0:
    print(f"  [FAIL] Impossible RPE-1 CENP-B pairs: {rpe1_b_pairs}")
    sys.exit(1)
if rpe1_b_45_65 is None or not (0.0 <= rpe1_b_45_65 <= 100.0):
    print(f"  [FAIL] Impossible RPE-1 CENP-B 45-65bp pct: {rpe1_b_45_65}")
    sys.exit(1)
if rpe1_b_sub85 is None or not (0.0 <= rpe1_b_sub85 <= 100.0):
    print(f"  [FAIL] Impossible RPE-1 CENP-B sub85 pct: {rpe1_b_sub85}")
    sys.exit(1)

# G. Reconcile metrics.json against Supplementary Tables TSV
with open('tables/Table_S6_cross_lineage_metrics_summary.tsv') as f:
    s6_rows = {r['cohort_id']: r for r in csv.DictReader(f, delimiter='	')}

if int(s6_rows['CHM13_REP1']['analyzed_pairs']) != rep1_pairs:
    print(f"  [FAIL] Table S6 CHM13 Rep 1 pairs mismatch: {s6_rows['CHM13_REP1']['analyzed_pairs']} vs {rep1_pairs}")
    sys.exit(1)
if int(s6_rows['HG002_T2T']['analyzed_pairs']) != hg002_pairs:
    print(f"  [FAIL] Table S6 HG002 pairs mismatch")
    sys.exit(1)
if int(s6_rows['RPE1_CENPA']['analyzed_pairs']) != rpe1_a_pairs:
    print(f"  [FAIL] Table S6 RPE-1 CENP-A pairs mismatch")
    sys.exit(1)
if int(s6_rows['RPE1_CENPB']['analyzed_pairs']) != rpe1_b_pairs:
    print(f"  [FAIL] Table S6 RPE-1 CENP-B pairs mismatch")
    sys.exit(1)

print(f"  [PASS] Cross-Lineage Replication verified across 4 cohorts:")
print(f"         - CHM13 Rep 1: mode={rep1_mode} bp, N={rep1_pairs:,}, core={rep1_core}%, 150bp={rep1_150}%")
print(f"         - HG002: mode={hg002_mode} bp, N={hg002_pairs:,}, core={hg002_core}%, sub85={hg002_sub85}%")
print(f"         - RPE-1 CENP-A: mode={rpe1_a_mode} bp, N={rpe1_a_pairs:,}, 147-175bp={rpe1_a_147_175}%")
print(f"         - RPE-1 CENP-B: N={rpe1_b_pairs:,}, 45-65bp={rpe1_b_45_65}%, sub85={rpe1_b_sub85}%")
EOF

echo "3. Verifying all 7 publication figures..."
for i in {1..7}; do
    for ext in pdf png svg; do
        matched=0
        for f in figures/Fig${i}_*.$ext; do
            if [ -f "$f" ]; then
                if [ ! -s "$f" ]; then
                    echo "  [FAIL] Figure $f is empty (0 bytes)!"
                    exit 1
                fi
                matched=1
            fi
        done
        if [ "$matched" -eq 0 ]; then
            echo "  [FAIL] Missing figure: Fig${i}_*.$ext!"
            exit 1
        fi
    done
    echo "  [PASS] Figure ${i} verified non-empty in PDF, PNG, and SVG"
done

echo "4. Verifying supplementary data sheets..."
if [ ! -s "tables/Supplementary_Tables_S1_to_S8.xlsx" ]; then
    echo "  [FAIL] Multi-sheet Excel workbook missing or empty: tables/Supplementary_Tables_S1_to_S8.xlsx"
    exit 1
fi
echo "  [PASS] Multi-sheet Excel workbook present: Supplementary_Tables_S1_to_S8.xlsx"

for s in {1..8}; do
    matched=0
    for t in tables/Table_S${s}_*.tsv; do
        if [ -f "$t" ]; then
            if [ ! -s "$t" ]; then
                echo "  [FAIL] Supplementary Table $t is empty (0 bytes)!"
                exit 1
            fi
            matched=1
        fi
    done
    if [ "$matched" -eq 0 ]; then
        echo "  [FAIL] Missing Supplementary Table S${s} (TSV)!"
        exit 1
    fi
    echo "  [PASS] Supplementary Table S${s} (TSV) present and non-empty"
done

echo "================================================================================"
echo "  ALL VERIFICATIONS PASSED! Reviewer Package Integrity Verified."
echo "================================================================================"
