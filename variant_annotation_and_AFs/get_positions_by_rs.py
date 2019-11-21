from cruzdb import Genome
import sys

"""
Gets SNP chr:position by rs id
"""

fname = sys.argv[1]
snp_col = sys.argv[2]
hg38 = Genome('hg38')
dbsnp = hg38.snp151

for rs in (l.split()[0].rstrip() for l in open(fname) if l.startswith("rs")):
with open(fname) as f:
	for l in f:
		spl = l.strip().split("\t")
		rs = spl[snp_col]
	    snp = dbsnp.filter_by(name=rs).first()
    	if snp:
        	spl[snp_col] = str(snp.chromStart) + ":" + str(snp.chromEnd)
    	print "\t".join(spl)
