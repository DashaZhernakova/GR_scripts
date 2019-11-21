#!/bin/bash

#
# Runs a liftover/liftdown for a plink genotype dataset
# run as: listover.sh <plink_bfile_base> <abbreviation of the target genome build> <liftover chain file>
#
# Requires UCSC liftover tool available at: http://hgdownload.soe.ucsc.edu/admin/exe/

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

x=`wc -l liftover_unmapped`
echo
echo "`echo $x | awk '{print $1/2}'` variants failed the liftover"
echo

#
# Update positions in the plink bim file after liftover
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

plink \
--bfile ${plink_f}.${to}.sorted \
--missing \
--allow-extra-chr \
--out ${plink_f}.${to}.sorted

#
# Remove variants with same positions (keep the one with the lowest missing rate)
#
~/scripts/exclude_duplicate_pos_in_bim_missingness.py  ${plink_f}.${to}.sorted.lmiss ${plink_f}.${to}.sorted.bim > ${plink_f}.${to}.sorted.noduplicates

nodup=`echo $(wc -l ${plink_f}.${to}.sorted.noduplicates) | awk '{print $1}'`
novar=`echo $(wc -l ${plink_f}.${to}.sorted.bim) | awk '{print $1}'`
echo
echo "Found $(($novar-$nodup)) duplicates"
echo

plink \
--make-bed \
--bfile ${plink_f}.${to}.sorted \
--allow-extra-chr \
--extract ${plink_f}.${to}.sorted.noduplicates \
--out ${plink_f}.${to}.sorted.rmdup

rm ${plink_f}.positions.bed
rm snps.txt
rm ${plink_f}.${to}.bim.tmp
rm ${plink_f}.${to}.bim
rm ${plink_f}.${to}.fam
rm ${plink_f}.${to}.sorted.bed
rm ${plink_f}.${to}.sorted.bim
rm ${plink_f}.${to}.sorted.fam
rm ${plink_f}.${to}.sorted.duplicates
rm ${plink_f}.${to}.sorted.*miss
rm *.nosex