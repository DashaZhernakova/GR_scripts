#!/usr/local/bin/python
import sys

iid_dict = {}
fid_dict = {}
with open(sys.argv[2]) as sample_f:
    for l in sample_f:
        spl = l.strip().split("\t")
        iid_dict[spl[1]] = spl[3]
        fid_dict[spl[0]] = spl[2]
        
with open(sys.argv[1]) as fam_f:
    for l in fam_f:
        spl = l.strip().split()
        
        iid = iid_dict[spl[1]]
        fid = fid_dict[spl[0]]
        pid = "0"
        mid = "0"

        if spl[2] != "0":
            pid = iid_dict[spl[2]]
        if spl[3] != "0":
            mid = iid_dict[spl[3]]
            
        print "\t".join([fid, iid, pid, mid, spl[4], spl[5]])