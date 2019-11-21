#!/usr/local/bin/python

"""
Rename samples from original sample ids to "Source_Population_index"
"""
from collections import defaultdict
out = open("/Users/dashazhernakova/Documents/Doby/GenomeRussia/freeze2019/population_genetics/GR+papers_samples_2.txt", "w")

source_dict = {"Pagani" : "P", "Mallick" : "M", "GR" : "G"}
indices = defaultdict(int)
with open("/Users/dashazhernakova/Documents/Doby/GenomeRussia/freeze2019/population_genetics/GR+papers_samples.txt") as f:
    out.write(f.readline().strip() + "\tnew_id\n")
    for l in f:
        spl = l.strip().split("\t")
        pop = spl[2]
        if pop not in indices:
            indices[pop] = 1
        new_id = source_dict[spl[6]] + "_" + spl[2].replace("-", "").replace("_", "") + "_" + str(indices[pop])
        out.write("\t".join(spl) + "\t" + new_id + "\n")
        indices[pop] += 1

out.close()
        
