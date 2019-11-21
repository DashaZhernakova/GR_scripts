#!/usr/bin/python
import re
import sys

with open(sys.argv[1]) as f:
    print f.readline().strip()
    for l in f:
        spl = l.strip().split("\t")
        omim = "NA"
        if "OMIM" in spl[5]:
            m = re.search("OMIM_([0-9]+)_?", spl[5])
            omim = m.group(1)
        print "\t".join(spl[:6]) + "\t" + omim + "\t" + "\t".join(spl[6:])