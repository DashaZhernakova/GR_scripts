#!/usr/bin/python
from collections import defaultdict
import sys
import vcf

def addSNPsToDict(fname,snp_dict, idx):
        vcf_reader = vcf.Reader(open(fname))
        for record in vcf_reader:
                id = record.CHROM + ":" + str(record.POS)
                if 'MAF' not in record.INFO:
			print "no MAF", id
			continue
		maf = record.INFO['MAF'][0]
                ac = int(record.INFO['AC_Het'][0]) + int(record.INFO['AC_Hom'][0])
                an = int(record.INFO['NS'])
                snp = snp_dict[id]
                snp[idx] = maf
                snp[3] += ac
                snp[4] += an

        return snp_dict

def calcOverlaps(af_list,res):
	#p,n,y,pn,ny,py,pny
	if af_list[0] > 0:
		res[0] += 1
	if af_list[1] > 0:
                res[1] += 1
	if af_list[2] > 0:
                res[2] += 1
	if (af_list[0] > 0) and ( af_list[1] > 0):
		res[3] += 1
	if (af_list[1] > 0) and ( af_list[2] > 0):
                res[4] += 1
	if (af_list[0] > 0) and ( af_list[2] > 0):
                res[5] += 1
	if (af_list[0] > 0) and ( af_list[1] > 0) and  ( af_list[2] > 0):
		res[6] += 1
	return res


snp_dict = defaultdict(lambda: [0]*5)
snp_dict = addSNPsToDict(sys.argv[1],snp_dict, 0)
print "first done"
snp_dict = addSNPsToDict(sys.argv[2],snp_dict, 1)
print "second done"
snp_dict = addSNPsToDict(sys.argv[3],snp_dict, 2)
print "third done"

out = open(sys.argv[4],"w")
res = [0,0,0,0,0,0,0]
for snp, af_list in snp_dict.items():
        maf = af_list[3]/float(2*af_list[4])
	out.write(snp + "\t" + "\t".join(map(str,af_list[0:3])) + "\t" + str(maf) + "\n")
	res = calcOverlaps(af_list, res)
out.close()
print res
