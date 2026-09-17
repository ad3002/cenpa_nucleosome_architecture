#!/usr/bin/env python3
"""
sync_paired.py

Synchronizes paired-end FASTQ.gz files record-by-record, verifying header IDs.
Supports unlimited reads when max_reads <= 0.
Fails with non-zero exit code if input is empty, truncated, or headers do not match.
"""

import gzip
import sys
import os

if len(sys.argv) < 6:
    print("Usage: sync_paired.py <f1_in> <f2_in> <f1_out> <f2_out> <max_reads>")
    sys.exit(1)

f1_in, f2_in, f1_out, f2_out, max_reads_str = sys.argv[1:6]
max_reads = int(max_reads_str)

if not os.path.exists(f1_in):
    raise FileNotFoundError(f"Input file not found: {f1_in}")
if not os.path.exists(f2_in):
    raise FileNotFoundError(f"Input file not found: {f2_in}")

limit_msg = f"up to {max_reads:,} pairs" if max_reads > 0 else "unlimited pairs"
print(f"Syncing {f1_in} and {f2_in} ({limit_msg})...")

def open_gz_safe(p):
    return gzip.open(p, 'rt', errors='replace')

count = 0
os.makedirs(os.path.dirname(os.path.abspath(f1_out)), exist_ok=True)
os.makedirs(os.path.dirname(os.path.abspath(f2_out)), exist_ok=True)

with open_gz_safe(f1_in) as r1, open_gz_safe(f2_in) as r2, \
     gzip.open(f1_out, 'wt') as w1, gzip.open(f2_out, 'wt') as w2:
    while max_reads <= 0 or count < max_reads:
        h1 = r1.readline()
        h2 = r2.readline()
        
        # End of stream
        if not h1 and not h2:
            break
        if not h1 or not h2:
            raise ValueError(f"Premature truncation: one file ended before the other at record {count}.")

        s1 = r1.readline()
        p1 = r1.readline()
        q1 = r1.readline()

        s2 = r2.readline()
        p2 = r2.readline()
        q2 = r2.readline()

        if not (s1 and p1 and q1 and s2 and p2 and q2):
            raise ValueError(f"Truncated FASTQ record at record {count}.")

        id1 = h1.split()[0]
        id2 = h2.split()[0]
        # Clean common /1 and /2 suffixes if present
        b1 = id1[:-2] if id1.endswith(('/1', '.1')) else id1
        b2 = id2[:-2] if id2.endswith(('/2', '.2')) else id2

        if b1 != b2:
            raise ValueError(f"Mismatched paired-end ID at record {count}: {id1} vs {id2}")

        w1.write(h1 + s1 + p1 + q1)
        w2.write(h2 + s2 + p2 + q2)
        count += 1
        if count % 1000000 == 0:
            print(f"  Synced {count:,} pairs...")

if count == 0:
    raise ValueError(f"Zero valid paired records found in {f1_in} and {f2_in}.")

print(f"Successfully wrote {count:,} synchronized read pairs to {f1_out} and {f2_out}")
