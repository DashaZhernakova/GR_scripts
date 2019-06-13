#!/bin/bash

plink_f=$1
to=$2
chain_f=$3
awk 'BEGIN {FS=" "}; {OFS="\t"}; {print $1, $4-1, $4, $2}' \
${plink_f}.bim > ${plink_f}.positions.bed

#
# Liftover
#
~/tools/liftOver \
${plink_f}.positions.bed \
${chain_f} \
${plink_f}.positions.${to}.bed \
liftover_unmapped

#
# Update positions after liftover
#
cut -f4 ${plink_f}.positions.${to}.bed > snps.txt
plink \
--bfile ${plink_f} \
--extract snps.txt --make-bed \
--out ${plink_f}.${to}

mv ${plink_f}.${to}.bim ${plink_f}.${to}.bim.tmp
python ~/scripts/update_pos_after_liftover.py \
${plink_f}.positions.${to}.bed \
${plink_f}.${to}.bim.tmp \
> ${plink_f}.${to}.bim

rm ${plink_f}.positions.${to}.bed
#
# Sort the resulting file
#
plink \
--bfile ${plink_f}.${to} \
--make-bed --out ${plink_f}.${to}.sorted \
--allow-extra-chr

#
# Remove variants with same positions (keep the 1st)
#
plink \
--bfile ${plink_f}.${to}.sorted \
--allow-extra-chr \
--list-duplicate-vars suppress-first

plink \
--make-bed \
--bfile ${plink_f}.${to}.sorted \
--allow-extra-chr \
--exclude plink.dupvar \
--out ${plink_f}.${to}.sorted.rmdup