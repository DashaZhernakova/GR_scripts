#!/usr/bin/python
import sys
import vcf
from geneimpacts import VEP, Effect
from collections import defaultdict

def getTopImpact(effects):
	top_severe = Effect.top_severity(effects)
	if type(top_severe) is list:
		top_impact = top_severe[0].top_consequence
	else:
		top_impact = top_severe.top_consequence
	return top_impact

def getGenesAndTranscrPerImpact(effects):
	gene_set = set()
	impact_dict = defaultdict(int)
	novel = '1'
	for effect in effects:
		gene = effect.effects['SYMBOL']
		if len(gene) > 1:
			gene_set.add(effect.effects['SYMBOL'])
		impact_dict[effect.impact_severity] += 1
		if effect.effects['Existing_variation']:
			novel = '0'
	if len(gene_set) == 0:
		gene_set.add('intergenic')
	return gene_set,impact_dict, novel

def isLoF(effects):
	for effect in effects:
		if effect.effects['LoF'] == 'HC':
			return '1'
	return '0'

def getGenotypeCounts(record, samples):
	samples_gen_counts = [0, 0, 0, 0, 0]

	for s in samples:
		#print s
		gen = record.genotype(s).gt_type
		#print s, gen
		if gen >= 0:
			samples_gen_counts[gen] += 1
			print samples_gen_counts
		else:
			samples_gen_counts[3] += 1
	sum_counts = 2*sum(samples_gen_counts[0:3])
	if sum_counts > 0:
		samples_gen_counts[4] = float(samples_gen_counts[1] + 2*samples_gen_counts[2])/sum_counts
	return samples_gen_counts

def loadCustomAnnotation(fname, col1, col2):
	f = open(fname)
	annot_dict = {}
	for l in f:
		spl = l.rstrip().split("\t")
		annot_dict[spl[col1]] = spl[col2]
	return annot_dict

def lookupGeneInAnnotation(gene_set, annot):
	annot_str = ""
	for g in gene_set:
		if g in annot:
			annot_str += ";" + annot[g]
	return annot_str.replace(";", "", 1)

def lookup_CSQ(effects, field):
	res_set = set()
	for effect in effects:
		if effect.effects[field] != "":
			res_set.add(effect.effects[field])
	if len(res_set) == 0:
		return "NA"
	return ";".join(res_set)

def getGroupSampleNames(fname, factor = False, col1 = 0, col2 = 0):
	samples = []
	with open(fname) as f:
		for l in f:
			spl = l.strip().split("\t")
			if factor:
				if spl[col2] == factor:
					samples.append(spl[col1])
				if factor == "all":
					samples.append(spl[col1])
			else:
				samples.append(spl[col1])
	return samples
def fillSampleMAFs(group_names, record):
	sample_descr_file = "/home/dzhernakova/GR/insert_size_800/vcf/ps+nov+yak/GR_samples.txt" 
	sample_mafs = {}
	for gr in group_names:
		print gr
		samples = getGroupSampleNames(sample_descr_file, gr, 0, 1)
		samples_gen_counts = getGenotypeCounts(record, samples)
		sample_mafs[gr] = samples_gen_counts[4]
		#print samples
		#print samples_gen_counts
	return sample_mafs

if __name__ == '__main__':
	vcf_fname = sys.argv[1]
	#vcf_fname = "/Users/dashazhernakova/Documents/Doby/GenomeRussia/60samples/ps+nov+yak+papers.filtered.from_russia.rsids.VEP.coding.vcf.gz"
	vcf_reader = vcf.Reader(open(vcf_fname))
	VEP.keys = vcf_reader.infos['CSQ'].desc.split(" ")[-1].split("|")
	
	sample_groups = ['Pskov', 'Novgorod', 'Yakut', 'all']
	
	
	trait = "NA"
	hgmd_annot = "NA"
	print "\t".join(["CHROM", "POS", "ID", "REF", "ALT", "genes", 
	"num_transcripts", "num_transcripts_HIGH_impact","num_transcripts_MED_impact","num_transcripts_LOW_impact",
	"top_impact"]) + "\t" + "\t".join(sample_groups) + "\t"+ "\t".join(["is LoF", "is novel", 
	"var associated trait", "gene associated trait", "HGMD_variant_annotation", 
	"1000G_EUR_MAF", "ExAC_EUR_MAF"])

	gwas_catalog = loadCustomAnnotation("/home/dzhernakova/resources/gwas_catalog_v1.0.1_e86_r2016-10-23.cut.gene_to_traits.txt", 0, 1)
	#gwas_catalog = loadCustomAnnotation("/Users/dashazhernakova/Documents/resources/hg38/GWAS_catalog/gwas_catalog_v1.0.1_e86_r2016-10-23.cut.gene_to_traits.txt", 0, 1)
	
	for record in vcf_reader:
		if not record.ID:
			record.ID = record.CHROM + ":" + str(record.POS)
		effects = []
		for annot in record.INFO['CSQ']:
			vep = VEP(annot)
			effects.append(vep)
		gene_set,impact_dict, novel = getGenesAndTranscrPerImpact(effects)
		sample_mafs = fillSampleMAFs(sample_groups, record)
		top_impact = getTopImpact(effects)
		
		
		kg_maf = lookup_CSQ(effects, "EUR_MAF")
		exac_maf = lookup_CSQ(effects, "ExAC_NFE_MAF")
	   
		trait = ""
		hgmd_annot = "" 
		if 'TRAIT' in record.INFO:
			trait = record.INFO['TRAIT']
		is_lof = isLoF(effects)
		
		if 'PHEN' in record.INFO:
			hgmd_annot = record.INFO['PHEN']
		
		
		maf_line = "\t".join([str(sample_mafs[gr]) for gr in sample_groups])
		print "\t".join([record.CHROM, str(record.POS), record.ID, record.REF, str(record.ALT), ",".join(gene_set),
		str(len(effects)), str(impact_dict["HIGH"]), str(impact_dict["MED"]), str(impact_dict["LOW"]),
		top_impact, maf_line, is_lof, novel, 
		trait, lookupGeneInAnnotation(gene_set, gwas_catalog), 
		hgmd_annot, kg_maf, exac_maf])
		break