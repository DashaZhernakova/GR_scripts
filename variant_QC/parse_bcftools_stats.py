#!/usr/bin/python
import os
import glob
import sys
from collections import defaultdict



def parse_stats_file(fname):
	"""
	Creates a dictionary with all bcftools stats from a stats file
	@param fname: output of bcftools stats command
	"""
    f = open(fname)
    bcftools_stats = dict()
    for line in f:
        s = line.split("\t")
        if s[0] == "ID":
            s_name = line.strip().split()[2].split("/")[0]            
            bcftools_stats = dict()
            bcftools_stats['novel'] = dict()
            bcftools_stats['known'] = dict()        

        if s[0] == "SN" and s[1] == '0':                    
            field = s[2].strip()[:-1]
            field = field.replace(' ', '_')
            value = float(s[3].strip())
            bcftools_stats['known'][field] = value
        if s[0] == "SN" and s[1] == '1':                   
            field = s[2].strip()[:-1]
            field = field.replace(' ', '_')
            value = float(s[3].strip())
            bcftools_stats['novel'][field] = value
        if s[0] == "TSTV" and s[1] == '0':
            fields = ['ts', 'tv', 'tstv']
            for i, f in enumerate(fields):
                value = float(s[i+2].strip())
                bcftools_stats['known'][f] = value
        if s[0] == "TSTV" and s[1] == '1':
            fields = ['ts', 'tv', 'tstv']
            for i, f in enumerate(fields):
                value = float(s[i+2].strip())
                bcftools_stats['novel'][f] = value
        if s[0] == "SiS" and s[1] == '0':
            fields = ['sis_allele_count', 'number_of_singletons']
            for i, f in enumerate(fields):
                value = float(s[i+2].strip())
                bcftools_stats['known'][f] = value
        if s[0] == "SiS" and s[1] == '1':
            fields = ['sis_allele_count', 'number_of_singletons']
            for i, f in enumerate(fields):
                value = float(s[i+2].strip())
                bcftools_stats['novel'][f] = value
    bcftools_stats['novel']['number_of_samples'] = bcftools_stats['known']['number_of_samples']
    return s_name, bcftools_stats

if __name__ == '__main__':
	"""
	Gathers QC statistics from per-population bcftools stats files and combined the QC results into one table
	"""
    path = sys.argv[1]
    bcftools_stats = defaultdict(list)
    fields = ['number_of_samples', 'number_of_records', 'number_of_SNPs', 'number_of_multiallelic_SNP_sites', 'number_of_singletons', 'tstv']
    for f in fields:
        bcftools_stats[f] = []
    population_list = []
    for infile in glob.glob( os.path.join(path, '*/*stats') ):
        sys.stderr.write("Extracting stats from " + infile + "\n")
        population, stats = parse_stats_file(infile)
        population_list.append(population)
             
        for f in fields:

            bcftools_stats[f].append(stats['known'][f])
            bcftools_stats[f].append(stats['novel'][f])
        
    print "\t" + "\t\t".join(population_list)
    print "\t" + "\t".join(['known\tnovel']*len(population_list))
    for f in fields:
        print f + "\t" + "\t".join(map(str,bcftools_stats[f]))
    
