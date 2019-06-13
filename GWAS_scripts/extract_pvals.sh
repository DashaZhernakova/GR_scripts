file_prefix=$1

dom_f=${file_prefix}.dominant.assoc.logistic.rsid.txt
rec_f=${file_prefix}.recessive.assoc.logistic.rsid.txt
allelic_f=${file_prefix}.allelic.assoc.logistic.rsid.txt
codom_f=${file_prefix}.codominant.assoc.logistic.rsid.txt


awk '{if ($5 == "TEST") print; if ($5 == "DOM") print}' $dom_f > ${dom_f}.pvals.txt
awk '{if ($5 == "TEST") print; if ($5 == "REC") print}' $rec_f > ${rec_f}.pvals.txt
awk '{if ($5 == "TEST") print; if ($5 == "ADD") print}' $allelic_f > ${allelic_f}.pvals.txt
awk '{if ($5 == "TEST") print; if ($5 == "ADD") print}' $codom_f > ${codom_f}.pvals.txt