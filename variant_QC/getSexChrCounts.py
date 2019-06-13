import vcf
import sys
from collections import defaultdict

if __name__ == '__main__':
    vcf_fname = sys.argv[1]
    vcf_reader = vcf.Reader(open(vcf_fname))
    chrX_dict_het  = defaultdict(int)
    chrX_dict_hom  = defaultdict(int)
    chrY_dict_het  = defaultdict(int)
    chrY_dict_hom  = defaultdict(int)

    for record in vcf_reader.fetch('X'):
        for s in record.samples:
            gt_type = s.gt_type
            if gt_type == 1:
                chrX_dict_het[s.sample] += 1
            elif gt_type == 0 or gt_type == 2:
                chrX_dict_hom[s.sample] += 1
    for record in vcf_reader.fetch('Y'):
        for s in record.samples:
            gt_type = s.gt_type
            if gt_type == 1:
                chrY_dict_het[s.sample] += 1
            elif gt_type == 0 or gt_type == 2:
                chrY_dict_hom[s.sample] += 1

    for s in chrX_dict_het:
        print s + "\t" + str(chrX_dict_het[s]) + "\t" + str(chrX_dict_hom[s]) + "\t" + str(chrY_dict_het[s]) + "\t" + str(chrY_dict_hom[s])
