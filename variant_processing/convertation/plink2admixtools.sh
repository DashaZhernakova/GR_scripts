#!/bin/bash

module load admixtools/5.1 

plink_f=$1
#sample_f="/mnt/genomerussia/dzhernakova/ancient_dna/samples_to_populations.txt"
sample_f=$2
pop_f="all_pops_quadro.txt"

scripts_dir=/mnt/genomerussia/projects/ancient_DNA/dzhernakova/admixtools_scripts/

awk '{$6="1"; print $0}' ${plink_f}.fam > ${plink_f}.fam.tmp
mv ${plink_f}.fam.tmp ${plink_f}.fam

sed "s:__PLINKFILE__:${plink_f}:g" ${scripts_dir}/param_convertf_plink_ancestrymap.TEMPLATE.txt | \
sed "s:__AMAPFILE__:${plink_f}:g" > param_convertf_plink_ancestrymap.txt

convertf -p param_convertf_plink_ancestrymap.txt

python ${scripts_dir}/update_ind_file_with_populations.py ${plink_f}.ind $sample_f > ${plink_f}_2.ind

sed "s:__AMAPFILE__:${plink_f}:g" ${scripts_dir}/param_dstats.TEMPLATE.txt | \
sed "s:__POPFILE__:${pop_f}:g" > param_dstats.txt

#/mnt/genomerussia/dzhernakova/ancient_dna/admixtools_scripts/make_all_possible_quadruples.py ${plink_f}.ind > all_pop_quadro.txt
