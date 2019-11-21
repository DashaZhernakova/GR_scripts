#!/usr/bin/python
import sys
import gzip

"""
Updates coordinates in bim file after liftover (part of the liftover.sh pipeline):
Writes new coordinates after liftover from a bed interval file resulting from liftover 
to the plink bim file on which the liftover was done.
SNP matching is done on SNP ids.
"""

bed = open(sys.argv[1]) # bed interval file with coordinates lifted over
bim = open(sys.argv[2]) # bim file to update

alleles = [l.strip().split("\t") for l in bim.readlines()]
bim.close()

line_num = 0
for l in bed:
  spl = l.strip().split("\t")
  bim_line = alleles[line_num]
  snp1 = bim_line[1]
  snp2 = spl[3]
  if not snp1 == snp2:
    print "Different SNPs, exiting! ", snp1, snp2
    sys.exit(-1)
  bim_line[0] = spl[0]
  bim_line[3] = spl[2]
  print "\t".join(bim_line)
  line_num += 1
bed.close()

