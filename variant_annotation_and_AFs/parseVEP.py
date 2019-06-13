#!/usr/bin/python
import sys
import vcf
from geneimpacts import VEP, Effect
from collections import defaultdict

def isLoF(effects):
    for effect in effects:
        if effect.effects['LoF'] == 'HC':
            return '1'
    return '0'

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

def loadCustomAnnotation(fname, col1, col2, col3 = -9):
    f = open(fname)
    annot_dict = {}
    for l in f:
        spl = l.rstrip().split("\t")
        if col3 == -9:
            annot_dict[spl[col1]] = spl[col2]
        else:
            annot_dict[spl[col1] + ":" + spl[col2]] = spl[col3]
    return annot_dict

def lookupGeneInAnnotation(gene_set, annot):
    annot_str = ""
    for g in gene_set:
        if g in annot:
            annot_str += ";" + annot[g]
    if annot_str == "":
        return "NA"
    return annot_str.replace(";", "", 1)

def lookupVariantInAnnotation(chr_pos, annot):
    annot_str = ""
    if chr_pos in annot:
        annot_str += ";" + annot[chr_pos]
    if annot_str == "":
        return "NA"
    return annot_str.replace(";", "", 1)

def lookup_CSQ(effects, field):
    res_set = set()
    for effect in effects:
        if effect.effects[field] != "":
            res_set.add(effect.effects[field])
    if len(res_set) == 0:
        return "NA"
    return ";".join(res_set)

if __name__ == '__main__':
    vcf_fname = sys.argv[1]
    #vcf_fname = "/Users/dashazhernakova/Documents/Doby/GenomeRussia/60samples/annotation/ps+nov+yak+papers.filtered.from_russia.rsids.VEP.coding.vcf.gz"
    vcf_reader = vcf.Reader(open(vcf_fname))
    VEP.keys = vcf_reader.infos['CSQ'].desc.split(" ")[-1].split("|")
    
    #out_f = open("/Users/dashazhernakova/Documents/Doby/tmp.txt", "w")    
    out_f = open(sys.argv[2], "w")
    trait = "NA"
    hgmd_annot = "NA"

    #gwas_catalog = loadCustomAnnotation("/home/dzhernakova/resources/variant_annotation/gwas_catalog_31-07-18_cut_fmt_merged.txt", 0, 1, 3)
    #clinvar = loadCustomAnnotation("/Users/dashazhernakova/Documents/resources/hg38/variation/clinvar_snps_cut_31-07-18_merged.txt", 0, 1, 3)
    gwas_catalog = loadCustomAnnotation("/home/dzhernakova/resources/variant_annotation/gwas_catalog_31-07-18_cut_fmt_merged.txt", 0, 1, 3)
    clinvar = loadCustomAnnotation("/home/dzhernakova/resources/variant_annotation/clinvar_snps_cut_31-07-18_merged.txt", 0, 1, 3)
    cnt = 0
    out_f.write("\t".join(["chr", "pos", "snp_id", "snp_alleles","alt_allele", "gene", "gene_type", "impact_severity", "effect_consequence", "AA_change", "Polyphen_prediction", "Sift_prediction", "isLoF", "GWAS catalog", "ClinVar"]) + "\n")
    for record in vcf_reader:
        chrpos = record.CHROM + ":" + str(record.POS)
        
        if not record.ID:
            record.ID = chrpos
        effects = []
        for annot in record.INFO['CSQ']:
            vep = VEP(annot)
            effects.append(vep)
        gene_set,impact_dict, novel = getGenesAndTranscrPerImpact(effects)
        #top_impact = getTopImpact(effects)

        effects_sorted = sorted(effects, reverse = True)
        out_set = set()
        
        for effect in effects_sorted:
            gene = effect.effects['SYMBOL']
            imp = effect.impact_severity

            
            polyphen = "NA"
            if effect.effects["PolyPhen"] != None:
                polyphen = effect.effects["PolyPhen"]
            
            sift = "NA"
            if effect.effects["SIFT"] != None:
                sift = effect.effects["SIFT"]
            
            aa_change = effect.effects['HGVSp']
            if len(aa_change) == 0:
                aa_change = "NA"
            
            novel = 1
            if effect.effects['Existing_variation']:
                novel = '0'
            
            gene = effect.effects['SYMBOL']
            if len(gene_set) == 0:
                gene = 'intergenic'
            
            lof = "NA"
            if effect.effects['LoF'] == 'HC':
                lof = "LoF"

           
            if (effect.effects['Allele'], gene, effect.biotype, imp) not in out_set:
                out_f.write("\t".join([record.CHROM, str(record.POS), record.ID, record.REF + "/" + ",".join(map(str,record.ALT)), effect.effects['Allele'], gene, effect.biotype, imp, effect.consequence, aa_change,  polyphen, sift, lof, lookupVariantInAnnotation(chrpos, gwas_catalog), lookupVariantInAnnotation(chrpos, clinvar)]) + "\n")
            out_set.add((effect.effects['Allele'], gene, effect.biotype, imp))
        cnt += 1
        if cnt % 10000 == 0:
            print "Prosessed", cnt, "variants"
            
out_f.close()