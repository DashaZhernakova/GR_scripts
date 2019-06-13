#!/usr/bin/python
import sys

plink_f=sys.argv[1]
f2exclude = sys.argv[2]
snps2exclude = set()

with open(f2exclude) as f:
    for l in f:
        spl = l.strip().split()
        snps2exclude.add(spl[0] + ":" + spl[1])

with open(plink_f) as f:
    print f.readline().strip()
    for l in f:
        spl = l.strip().split()