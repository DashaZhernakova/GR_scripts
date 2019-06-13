#!/usr/bin/python
import sys
import gzip

bed = open(sys.argv[1])
bim = open(sys.argv[2])

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

