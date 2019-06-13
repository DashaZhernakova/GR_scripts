#!/bin/bash

file_prefix=$1

echo
echo "################" 
echo "Sort, subset, add genes to plink GWAS result files"
echo "################"
echo

if [ $# -eq 0 ]
  then
    echo "No arguments supplied"
    echo "Usage: sort_annotate_gwas_res.sh <input_file_prefix>"
    echo "Exiting"
    exit 1
fi

dom_f=${file_prefix}.dominant.assoc.logistic.rsid.txt
rec_f=${file_prefix}.recessive.assoc.logistic.rsid.txt
allelic_f=${file_prefix}.allelic.assoc.logistic.rsid.txt
codom_f=${file_prefix}.codominant.assoc.logistic.rsid.txt

echo
echo "File prefix: $file_prefix"
echo "Processing the following files:"
echo "Dominant: $dom_f"
echo "Recessive: $rec_f"
echo "Allelic: $allelic_f"
echo "Codominant: $codom_f"
echo
echo "Writing the output to the current directory"
echo


echo "Processing dominant"
awk '{if ($5 == "TEST") print; if ($5 == "DOM" && $9 < 0.05) print}' $dom_f | sort -k9,9g > ${dom_f}.p0.05.sorted.txt
awk 'BEGIN {OFS="\t"}; {if ($1 != "CHR") print $1, $3-1, $3, $2}' ${dom_f}.p0.05.sorted.txt | sort -k1,1 -k2,2n | closestBed -a - -b /home/dzhernakova/resources/b37/Ensembl_b37_v75.genes.sorted.bed.gz -d  | awk '{print $4, $8 ":" $9}' > ${dom_f}.snp_annot.txt
plink --annotate ${dom_f}.p0.05.sorted.txt attrib=${dom_f}.snp_annot.txt --out ${dom_f}.p0.05.sorted

echo "Processing recessive"
awk '{if ($5 == "TEST") print; if ($5 == "REC" && $9 < 0.05) print}' $rec_f | sort -k9,9g > ${rec_f}.p0.05.sorted.txt
awk 'BEGIN {OFS="\t"}; {if ($1 != "CHR") print $1, $3-1, $3, $2}' ${rec_f}.p0.05.sorted.txt | sort -k1,1 -k2,2n | closestBed -a - -b /home/dzhernakova/resources/b37/Ensembl_b37_v75.genes.sorted.bed.gz -d  | awk '{print $4, $8 ":" $9}' > ${rec_f}.snp_annot.txt
plink --annotate ${rec_f}.p0.05.sorted.txt attrib=${rec_f}.snp_annot.txt --out ${rec_f}.p0.05.sorted

echo "Processing allelic"
awk '{if ($5 == "TEST") print; if ($5 == "ADD" && $9 < 0.05) print}' $allelic_f | sort -k9,9g > ${allelic_f}.p0.05.sorted.txt
awk 'BEGIN {OFS="\t"}; {if ($1 != "CHR") print $1, $3-1, $3, $2}' ${allelic_f}.p0.05.sorted.txt | sort -k1,1 -k2,2n | closestBed -a - -b /home/dzhernakova/resources/b37/Ensembl_b37_v75.genes.sorted.bed.gz -d  | awk '{print $4, $8 ":" $9}' > ${allelic_f}.snp_annot.txt
plink --annotate ${allelic_f}.p0.05.sorted.txt attrib=${allelic_f}.snp_annot.txt --out ${allelic_f}.p0.05.sorted

echo "Processing codominant"
awk '{if ($5 == "TEST") print; if ($5 == "ADD" && $9 < 0.05) print}' $codom_f | sort -k9,9g > ${codom_f}.p0.05.sorted.txt
awk 'BEGIN {OFS="\t"}; {if ($1 != "CHR") print $1, $3-1, $3, $2}' ${codom_f}.p0.05.sorted.txt | sort -k1,1 -k2,2n | closestBed -a - -b /home/dzhernakova/resources/b37/Ensembl_b37_v75.genes.sorted.bed.gz -d  | awk '{print $4, $8 ":" $9}' > ${codom_f}.snp_annot.txt
plink --annotate ${codom_f}.p0.05.sorted.txt attrib=${codom_f}.snp_annot.txt --out ${codom_f}.p0.05.sorted

echo
echo "Combining files into one: ${file_prefix}.all_tests.txt"
echo "NB! Codominant results are NOT included to avoid super-low p-values"
head -1 ${allelic_f}.p0.05.sorted.annot > ${file_prefix}.all_tests.txt
tail -n+2 ${allelic_f}.p0.05.sorted.annot | awk 'BEGIN {OFS="\t"}; {$5="ALLELIC"; print $0}' >> ${file_prefix}.all_tests.txt
tail -n+2 ${dom_f}.p0.05.sorted.annot | awk 'BEGIN {OFS="\t"}; {$5="DOM"; print $0}' >> ${file_prefix}.all_tests.txt
tail -n+2 ${rec_f}.p0.05.sorted.annot | awk 'BEGIN {OFS="\t"}; {$5="REC"; print $0}' >> ${file_prefix}.all_tests.txt
sort -k9,9g ${file_prefix}.all_tests.txt | uniq > ${file_prefix}.all_tests.sorted.txt

echo "Finished!"
