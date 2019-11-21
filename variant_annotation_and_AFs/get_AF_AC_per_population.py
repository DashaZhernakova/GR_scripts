#!/usr/bin/python
import sys
import vcf
from vcf.parser import _Info as VcfInfo, field_counts as vcf_field_counts
from collections import defaultdict
from collections import namedtuple

"""
Counts allele number and allele counts for each SNP in each population group.
"""

def getGroupSampleNames(fname, col1 = 0, col2 = 1):
    """
    Makes a dictionalry with population : list of samples in this population
    @param fname: name of the file containing sample id in col1 and population name in col2
    """
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
    @param record: VCF reader record to get counts from
    @param sample_popgroups: dictionary with population : list of samples in this population
    
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

def addMAFs2INFO(record, per_allele_counts):
	""" Adds AFs and ACs per population to the INFO field of the VCF file """
    infos = defaultdict(list)
    for per_pop_all_counts in per_allele_counts:
        for pop, count_list in per_pop_all_counts.items():
            infos[pop].append(count_list[2])
    for pop, pop_AF_list in infos.items():
        record.add_info('AF_' + pop, pop_AF_list)
    return record

def add_info(reader, info_records):
    """Add INFO records to a pyVCF reader."""
    for k in info_records:
        reader.infos[k.id] = k
    return reader
def new_info(tag, num, val_type, desc, source, version):
    """Create a VCF INFO record suitable for a pyVCF reader."""
    Info = namedtuple(
        'Info',
        ['id', 'num', 'type', 'desc', 'source', 'version'])
    return Info(tag, num, val_type, desc, source, version)

if __name__ == '__main__':    
    """
	Counts allele number and allele counts for each SNP in each population group.
	Outputs this as a tab-delimited table and optionally (if the third command line argument is given) 
	as a new VCF file with ACs, AFs and ANs per population written in the INFO field.
	sys.argv[1]: input vcf fname with genotypes
	sys.argv[2]: sample annotation file: sample \t population (if other format, change the column numbers)
	sys.argv[3] (optional): new VCF file name to write the VCF with AFs
	"""
    
    vcf_fname = sys.argv[1]
    vcf_reader = vcf.Reader(open(vcf_fname))
         
    sample_fname = sys.argv[2]
    sample_popgroups = getGroupSampleNames(sample_fname)
    population_lst = sample_popgroups.keys()
    header = "\t".join(["chr", "pos", "id", "ref", "alt", "gene", "trait"])
    new_infos = []
    for pop in population_lst:
        header += "\tAN_" + pop + "\tAC_"  + pop + "\tAF_"  + pop + "\tnHomo_"  + pop
        new_infos.append(new_info("AF_"  + pop, vcf_field_counts['A'], 'Float', 'Alternative allele frequency for ' + pop, None, None))
    
    write_new_vcf = False
    if len(sys.argv) > 3:
    	write_new_vcf = True
    	out_vcf_fname = sys.argv[3]
    	vcf_writer = vcf.Writer(open(out_vcf_fname, 'w'), add_info(vcf_reader, new_infos))
    	#vcf_writer = vcf.Writer(sys.stdout, add_info(vcf_reader, new_infos))
    
    
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
        if write_new_vcf:
	        new_record = addMAFs2INFO(record, per_allele_counts)
    	    vcf_writer.write_record(new_record)
       
    if write_new_vcf:
	    vcf_writer.close()