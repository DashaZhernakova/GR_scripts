#!/bin/bash

# Run ibd calculations
# input phased file: ../GR+papers.eurasian.chr${chr}.geno0.05.hwe1e-4.me0.05.maf0.01.phased.vcf.gz
# samples2exclude.txt contains children

for chr in `seq 1 22`
do
 #chr=$1
 echo "chr=$chr"
 java -Xmx100g -jar ~/tools/refined-ibd.12Jul18.a0b.jar \
 gt=../GR+papers.eurasian.chr${chr}.geno0.05.hwe1e-4.me0.05.maf0.01.phased.vcf.gz \
 out=chr${chr} \
 map=/home/dzhernakova/resources/recombination_maps/plink.chr${chr}.GRCh38.map \
 excludesamples=./samples2exclude.txt
done

rm -f GR+papers.eurasian.filtered.ibd
for chr in `seq 1 22`; do zcat chr${chr}.ibd.gz >> GR+papers.eurasian.filtered.ibd; done

# prepare plink genotype files for the subsequent script
#bcftools concat -Oz -o ../GR+papers.eurasian.all_chr.geno0.05.hwe1e-4.me0.05.maf0.01.phased.vcf.gz ../GR+papers.eurasian.chr*.geno0.05.hwe1e-4.me0.05.maf0.01.phased.vcf.gz
plink --vcf ../GR+papers.eurasian.all_chr.geno0.05.hwe1e-4.me0.05.maf0.01.phased.vcf.gz --make-bed --out GR+papers.eurasian.filtered --double-id

#prepare sample annotation 
cut -d " " -f 2 GR+papers.eurasian.filtered.fam > sample_list.txt
~/scripts/grepf.py -f /mnt/sochi/dzhernakova/GR/GR_freeze2019/genotypes/merged/GR+papers/sample2pop_renamed.txt -f_col 0 -i sample_list.txt -i_col 0 --f_add | cut -f2,3 > sample_annotation.txt

c=`wc -l sample_list.txt`
x=`echo $c | awk '{print $1}'`

c=`wc -l sample_annotation.txt`
y=`echo $c | awk '{print $1}'`


if [[ "$x" != "$y" ]]; then
 echo "ERROR: Something went wrong with sample annotation - check the number of samples in the fam file and in the annotation"
 echo "Exiting."
 exit 1
fi

# aggregate ibd results together, normalize by population size. Script by Lena Lukianova
# output file .segmentcounts - overal number of shared segments
# output file .segments - number of shared segments normalized by population size. This file should be used for plotting and other analyses!
perl ~/scripts/crimean_ibd.pl dataset_fname=GR+papers.eurasian.filtered.ibd anno_fname=sample_annotation.txt snps_ind=0 fam_name=GR+papers.eurasian.filtered

sed -i 's:,:.:g' GR+papers.eurasian.filtered.ibd.format.segments.step0.5