#/usr/bin/python
import sys
from collections import defaultdict
import gzip


def getCallRate(mis_fname, bim_fname):
	"""
	Finds duplicate SNPs in a plink genotype fileset and prints SNPs with higher missing rate to be excluded
	(as part of the liftover.sh pipeline)
	@param fname: file containing missing rate per SNP (.lmiss file, output of plink --missing command)
	@param bim_f: plink bim file name, on which missingness was calulated
	"""

	mis_f = open(mis_fname)
	bim_f = open(bim_fname)
	cur_mis = 1
	cur_id = ""
	cur_spl = ["","","","","",""]
	cnt = 0
	mis_f.readline()
	for l in bim_f.readlines():
		bim_spl = l.strip().split()
		mis_spl= mis_f.readline().strip().split()
		if cur_spl == bim_spl:
			cur_mis = 1
            cur_id = ""
            cur_spl = ["","","","","",""]
			continue
		if bim_spl[0] + ":" + bim_spl[3] != cur_spl[0] + ":" + cur_spl[3]: # new SNP
			if cur_id != "": print cur_id
			cur_mis = float(mis_spl[4])
			cur_id = bim_spl[1]
			cur_spl = bim_spl
		elif float(mis_spl[4]) <= cur_mis: # same with lower missingness
			cur_mis = float(mis_spl[4])
			cur_id = bim_spl[1]
			cur_spl = bim_spl

			
getCallRate(sys.argv[1], sys.argv[2])