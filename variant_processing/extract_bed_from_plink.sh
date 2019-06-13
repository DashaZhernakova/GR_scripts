#!/bin/bash

in_plink=$1
mask=$2
out_plink=$3

echo "Filtering $mask from $in_plink"

awk 'BEGIN {FS="\t"}; {OFS="\t"}; {print $1, $4-1, $4, $2}' \
${in_plink}.bim \
> tmp.bed

intersectBed \
-a tmp.bed \
-b ${mask} \
| cut -f4 \
> tmp_filtered.txt

plink \
--bfile ${in_plink} \
--extract tmp_filtered.txt \
--make-bed \
--out ${out_plink}

rm tmp.bed
rm tmp_filtered.txt
echo "Finished"