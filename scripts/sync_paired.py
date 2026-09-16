#!/usr/bin/env python3
import gzip
import sys

f1_in, f2_in, f1_out, f2_out, max_reads = sys.argv[1:6]
max_reads = int(max_reads)

def open_gz_safe(p):
    return gzip.open(p, 'rt', errors='ignore')

print(f'Syncing {f1_in} and {f2_in} up to {max_reads} pairs...')

count = 0
with open_gz_safe(f1_in) as r1, open_gz_safe(f2_in) as r2, \
     gzip.open(f1_out, 'wt') as w1, gzip.open(f2_out, 'wt') as w2:
    while count < max_reads:
        try:
            h1 = r1.readline()
            s1 = r1.readline()
            p1 = r1.readline()
            q1 = r1.readline()
            
            h2 = r2.readline()
            s2 = r2.readline()
            p2 = r2.readline()
            q2 = r2.readline()
        except Exception:
            break
            
        if not (h1 and s1 and p1 and q1 and h2 and s2 and p2 and q2):
            break
            
        id1 = h1.split()[0]
        id2 = h2.split()[0]
        if id1 != id2:
            print(f'Mismatched pair at record {count}: {id1} vs {id2}')
            break
            
        w1.write(h1 + s1 + p1 + q1)
        w2.write(h2 + s2 + p2 + q2)
        count += 1
        if count % 1000000 == 0:
            print(f'  Synced {count:,} pairs...')

print(f'Successfully wrote {count:,} synchronized read pairs to {f1_out} and {f2_out}')
