#!/bin/bash

file_prefix=$1
out_folder=$2

echo
echo "################" 
echo "Add rs ids to plink GWAS results"
echo "################"
echo

if [ $# -eq 0 ]
  then
    echo "No arguments supplied"
    echo "Usage: add_rsids.sh <input_file_prefix>"
    echo "Exiting"
    exit 1
fi

echo "File prefix: $file_prefix"
echo "Adding rs ids to the following files: ${file_prefix}*assoc.logistic"
echo "Writing to the current directory"
echo

for f in ${file_prefix}*assoc.logistic
do
 fname=${f##*/}
 echo $fname
 /home/dzhernakova/scripts/add_rs_to_gwas_res.py $f /home/dzhernakova/resources/b37/variation/dbSNP_20160408.SNPs.vcf.gz > ${out_folder}/${fname}.rsid.txt
 returnCode=$?
 if [ ! $returnCode -eq 0 ]; then
	echo "Failed, exiting"
	exit 1
 fi
done

echo "Finished all"
