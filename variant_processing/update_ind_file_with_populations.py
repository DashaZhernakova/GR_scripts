#!/usr/bin/python
import sys

pops = {}
with open(sys.argv[2]) as pop_f:
	for l in pop_f:
		spl = l.strip().split("\t")
		pops[spl[0]] = spl[1]

with open(sys.argv[1]) as ind_f:
	for l in ind_f:
		spl = l.strip().split()
		if spl[0] not in pops:
			sys.stderr.write("No sample " + spl[0] + " in population file! Exiting!")
			sys.exit(1)
		pop = pops[spl[0]]
		print spl[0] + "\t" + spl[1] + "\t" + pop
