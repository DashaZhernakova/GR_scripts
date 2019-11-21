#!/usr/local/bin/python
from collections import defaultdict
fname = "/Users/dashazhernakova/Documents/Doby/GenomeRussia/freeze2019/HLA_typing/hla_subset_res.txt"
out = open("/Users/dashazhernakova/Documents/Doby/GenomeRussia/freeze2019/HLA_typing/hla_subset_res_AFs.txt", "w")

pop_ac_dict = defaultdict(int)
all_alleles = []
with open(fname) as f:
    for l in f:
        spl = l.strip().split("\t")
        population = spl[1]
        pop_ac_dict[population] = defaultdict(int)

with open(fname) as f:
    for l in f:
        spl = l.strip().split("\t")
        population = spl[1]
        for allele in spl[2:]:
            if allele == "?":
                continue
            gene = allele.split("*")[0]
            
            pop_ac_dict[population][gene] += 1
            pop_ac_dict[population][allele] += 1

            if allele not in all_alleles:
                all_alleles.append(allele)

all_alleles = sorted(all_alleles)
out.write("population\t" + "\t".join(all_alleles) + "\n")
for population, ac_dict in pop_ac_dict.items():   
    out_array = ["0"]*len(all_alleles)
    for allele, ac in ac_dict.items():
        gene = allele.split("*")[0]
        if ":" in allele:
            out_array[all_alleles.index(allele)] = str(float(ac)/ac_dict[gene])
            #out_array[all_alleles.index(allele)] = str(ac)

    out.write(population + "\t" + "\t".join(out_array) + "\n")