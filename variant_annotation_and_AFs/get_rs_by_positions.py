#!/usr/local/bin/python
import sys
import vcf 
from collections import defaultdict

def getRs(chr_pos, vcf_reader):
    chr = chr_pos.split(":")[0]
    pos = int(chr_pos.split(":")[1])
    for record in vcf_reader.fetch(chr, pos - 1, pos):
        #print "##", record.ID, record
        if record.is_snp:
            return record.ID
    return chr_pos

in_f = sys.argv[1]
vcf_reader = vcf.Reader(open(sys.argv[2]))
col_num = int(sys.argv[3])

snp2rs = defaultdict()
with open(in_f) as f:
    for l in f:
        spl = l.strip().split("\t")
        snp = spl[col_num]
        rs = snp2rs.get(snp, None)
        #print "#", snp, rs
        if not rs:
            rs = getRs(snp, vcf_reader)
            snp2rs[snp] = rs
        spl[col_num] = rs
        print "\t".join(spl)
