from collections import defaultdict
import vcf
from geneimpacts import VEP, Effect


def overlapsExisting(vep_effects):
    """
    Check whether the variant is found in 1000G populations. Doesn't yet take MAF into account.
    """
    overlap_eur = False
    overlap_afr = False
    overlap_amr = False
    overlap_eas = False
    overlap_sas = False
    for effect in vep_effects:
        if effect.effects['AFR_MAF']: overlap_afr = True
        if effect.effects['AMR_MAF']: overlap_amr = True
        if effect.effects['EAS_MAF']: overlap_eas = True
        if effect.effects['EUR_MAF']: overlap_eur = True
        if effect.effects['SAS_MAF']: overlap_sas = True

    return overlap_eur, overlap_afr, overlap_amr, overlap_eas, overlap_sas

def isLoF(vep_effects):
    for effect in vep_effects:
        if effect.effects['LoF'] == 'HC':
            return True

if __name__ == '__main__':
    vcf_fname = "/home/dzhernakova/test/GR/genotyping/vcf/batch1/batch1.sorted.minQUAL40.minGQ20.minDP10.maxSP20.rsids.VEP.vcf.gz"
    vcf_reader = vcf.Reader(open(vcf_fname))
    VEP.keys = vcf_reader.infos['CSQ'].desc.split(" ")[-1].split("|")

    novel = False
    overlap_effects = defaultdict(int)
    no_overlap_effects = defaultdict(int)

    for record in vcf_reader:
        maf = float(2 * record.num_hom_alt + record.num_het) / (2*record.num_called)
        effects = []
        for annot in record.INFO['CSQ']:
            vep = VEP(annot)
            effects.append(vep)
        effects = sorted(effects) # last record - the most severe

        top_severe = Effect.top_severity(effects)
        if type(top_severe) is list:
            top_impact = top_severe[0].impact_severity
        else:
            top_impact = top_severe.impact_severity
        is_lof = isLoF(effects)
        overlap_eur, overlap_afr, overlap_amr, overlap_eas, overlap_sas = overlapsExisting(effects)
        overlap_1kg = overlap_eur or overlap_afr or overlap_amr or overlap_eas or overlap_sas
        if not overlap_1kg:
            no_overlap_effects[top_impact] += 1
            if is_lof: no_overlap_effects['LoF'] += 1
        else:
            overlap_effects[top_impact] += 1
            if is_lof: overlap_effects['LoF'] += 1

        if overlap_eur: count_dict['eur'] += 1
        if overlap_afr: count_dict['afr'] += 1
        if overlap_amr: count_dict['amr'] += 1
        if overlap_eas: count_dict['eas'] += 1
        if overlap_sas: count_dict['sas'] += 1

    for k, v in count_dict.items():
        print k, v
    print "Overlaping:"
    for k,v in overlap_effects:
        print k, v
    print "NOT overlaping:"
    for k,v in no_overlap_effects:
        print k, v
