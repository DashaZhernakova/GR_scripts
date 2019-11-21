#!/usr/bin/python
import sys
from collections import defaultdict

"""
Creates a plink fam file from a tab-separated sample desciption table.
Table columns: Sample ID Region Location Ethnicity Family_ID Family_status Gender

!!!! BEWARE OF FAMILIES WITH 2 AND MORE CHILDREN!!! NOT PARSED CORRECTLY
"""


gender_dict = {}
family_dict = {}
with open(sys.argv[1]) as f:
    f.readline()
    for l in f:
        spl = l.strip().split("\t")
        id = spl[0]
        fid = spl[4]
        gender_dict[id] = spl[6]

        if fid == 'None' or fid == '':
            family_dict[id] = ['0','0',id]
            continue
        
        if fid not in family_dict:
            family_dict[fid] = ['0','0','0']
        if spl[5] == 'Father':
            family_dict[fid][0] = id
        elif spl[5] == 'Mother':
            family_dict[fid][1] = id
        elif spl[5] == 'Child':
            if family_dict[fid][2] == '0':
                family_dict[fid][2] = id
            else:    
                family_dict[fid][2] += "," + id # if >2 children, add child id to the same family
            
            



for fid,f_dict in family_dict.items():
    if f_dict[0] != '0':
        print "\t".join([fid, f_dict[0],'0', '0', gender_dict[f_dict[0]], "-9"])
    if f_dict[1] != '0':
        print "\t".join([fid, f_dict[1],'0', '0', gender_dict[f_dict[1]], "-9"])
    if f_dict[2] != '0':
        for ch_id in f_dict[2].split(","):
            print "\t".join([fid, ch_id,f_dict[0], f_dict[1], gender_dict[ch_id], "-9"])