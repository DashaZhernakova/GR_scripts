from cruzdb import Genome
import sys

fname = "/Users/dashazhernakova/Documents/Doby/GenomeRussia/freeze2019/pharmacogenetics/uniq_snp_pharm_list.txt"
out_fname = "/Users/dashazhernakova/Documents/Doby/GenomeRussia/freeze2019/pharmacogenetics/uniq_snp_pharm_list.bed"
out = open(out_fname, "w")
hg38 = Genome('hg38')
dbsnp = hg38.snp151

for rs in (l.split()[0].rstrip() for l in open(fname) if l.startswith("rs")):
    snp = dbsnp.filter_by(name=rs).first()
    if snp:
        out.write(snp.chrom + "\t" + str(snp.chromStart) + "\t" + str(snp.chromEnd) + "\t" + rs + "\n")
    else:
        print "No SNP:", rs

out.close()