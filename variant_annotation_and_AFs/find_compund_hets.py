#!/usr/bin/python
import sys
import vcf
from geneimpacts import VEP, Effect
from collections import defaultdict

def fillSampleMAFs(group_names, record):
	sample_descr_file = "/Users/dashazhernakova/Documents/Doby/GenomeRussia/60samples/simons+pagani_samples.txt"
	sample_mafs = {}
	for gr in group_names:
		samples = getGroupSampleNames(sample_descr_file, gr, 0, 1)
		samples_gen_counts = getGenotypeCounts(record, samples)
		sample_mafs[gr] = samples_gen_counts[4]
		
	return sample_mafs

def getGenotypeCounts(record, samples):
	samples_gen_counts = [0, 0, 0, 0, 0]
	for s in samples:
		gen = record.genotype(s).gt_type
		if gen:
			samples_gen_counts[gen] += 1
		else:
			samples_gen_counts[3] += 1
	samples_gen_counts[4] = float(samples_gen_counts[1] + 2*samples_gen_counts[2])/2*sum(samples_gen_counts[0:3])
	return samples_gen_counts

if __name__ == '__main__':
	#vcf_fname = sys.argv[1]
	vcf_fname = "/Users/dashazhernakova/Documents/Doby/GenomeRussia/60samples/ps+nov+yak+papers.filtered.from_russia.rsids.VEP.coding.vcf.gz"
	vcf_reader = vcf.Reader(open(vcf_fname))
	VEP.keys = vcf_reader.infos['CSQ'].desc.split(" ")[-1].split("|")
	
	sample_groups = ['Pskov', 'Novgorod', 'Yakut', 'Simons', 'Pagani', 'all']
	transcr_dict = defaultdict[list]
	
	for record in vcf_reader:
		if not record.ID:
			record.ID = record.CHROM + ":" + str(record.POS)
		effects = []
		for annot in record.INFO['CSQ']:
			vep = VEP(annot)
			
			if vep.effects['LoF'] == 'HC':
				sample_mafs = fillSampleMAFs(sample_groups, record)
				transcr_dict[vep.effects['Feature'] + ":" + vep.effects['Symbol']].append([record.ID, sample_mafs])
	
	for transcr, snps in transcr_dict.items():
		for snp in snps:
			print transcr.split(":")[1] + "\t" + transcr.split(":")[0] + "\t". join(snp)