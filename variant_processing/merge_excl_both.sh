#!/bin/bash

f1=$1
f2=$2
merged=$3

plink --bfile ${f1} --bmerge ${f2} --make-bed --out ${merged} --allow-no-sex
plink --bfile ${f2} --flip ${merged}-merge.missnp --make-bed --out ${f2}.flipped --allow-no-sex
plink --bfile ${f1} --bmerge ${f2}.flipped --make-bed --out ${merged} --allow-no-sex
plink --bfile ${f2}.flipped --exclude ${merged}-merge.missnp --make-bed --out ${f2}.flipped.snp_excl --allow-no-sex
plink --bfile ${f1} --exclude ${merged}-merge.missnp --make-bed --out ${f1}.snp_excl --allow-no-sex
plink --bfile ${f1}.snp_excl --bmerge ${f2}.flipped.snp_excl --make-bed --out ${merged} --allow-no-sex

rm ${f1}.snp_excl*
rm ${f2}.flipped*
rm ${merged}-merge.missnp

