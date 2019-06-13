#!/usr/bin/python
import sys
import vcf
from collections import defaultdict


def getGroupSampleNames(fname, col1 = 0, col2 = 1):
    sample_popgroups = defaultdict(list)
    with open(fname) as f:
        for l in f:
            spl = l.strip().split("\t")
            sample_popgroups[spl[col2]].append(spl[col1])
    return sample_popgroups

def getGenotypeCounts(record, sample_popgroups):
    """
    Fills allele number, alt allele count, alt allele frequency and number of homozygotes by alt allele in all individuals per population from samples dictionary.
    Returns a list of dictionaries, where the  list index corresponds to the alternative allele number; the dictionaries contain population mapped to a list containing quadruples [AN, AC, AF, nHomo] for each population
    """
    per_allele_counts = []
    #print record.ALT
    for alt_al_num in range(1,len(record.ALT) + 1):
        per_pop_counts = {}
        #print alt_al_num
        for pop, samples in sample_popgroups.items():
            
            samples_gen_counts = [0, 0, 0, 0]       #AN, AC, AF, nHomo for population pop
            for s in samples:
                
                gen = record.genotype(s)['GT']
                #print s, gen
                if not gen == './.':
                    samples_gen_counts[0] += 2      # overall number of called alleles
                if gen == str(alt_al_num) + "/" + str(alt_al_num):      # number of alternative homozygotes
                    samples_gen_counts[3] += 1
                samples_gen_counts[1] += gen.count(str(alt_al_num))      # number of current alternative alleles
            if samples_gen_counts[0]  > 0:
                samples_gen_counts[2] = float(samples_gen_counts[1])/samples_gen_counts[0]      # AF = AC/AN
            #print samples_gen_counts
            per_pop_counts[pop] = samples_gen_counts
        
        per_allele_counts.append(per_pop_counts)

    return per_allele_counts


if __name__ == '__main__':
    #sample_fname = "/Users/dashazhernakova/Documents/Doby/GenomeRussia/60samples/samples2group.txt"
    sample_fname = sys.argv[2]
    sample_popgroups = getGroupSampleNames(sample_fname)
    population_lst = sample_popgroups.keys()
    header = "\t".join(["chr", "pos", "id", "ref", "alt", "gene", "trait"])
    for pop in population_lst:
        header += "\tAN_" + pop + "\tAC_"  + pop + "\tAF_"  + pop + "\tnHomo_"  + pop
    print header
    cnt = 0
    #vcf_fname = "/Users/dashazhernakova/Documents/Doby/GenomeRussia/60samples/chr22.pskov+novgorod+yakutia.SNPs.filtered.vcf.gz"
    vcf_fname = sys.argv[1]
    vcf_reader = vcf.Reader(open(vcf_fname))
    for record in vcf_reader:
        if not record.ID:
			record.ID = record.CHROM + ":" + str(record.POS)
        trait = record.INFO.get("TRAIT", "NA")
        gene = record.INFO.get("GENE", "NA")
        per_allele_counts = getGenotypeCounts(record, sample_popgroups)
        for alt_al_num, alt_al in enumerate(record.ALT):
            per_pop_all_counts = per_allele_counts[alt_al_num]
            out_str = "\t".join([record.CHROM, str(record.POS), record.ID, str(record.REF), str(alt_al), gene, trait])
            for pop in population_lst:
                out_str += "\t" + "\t".join(map(str,per_pop_all_counts[pop]))
            print out_str
        cnt += 1
        #if cnt > 100:
        #    break        
    
    