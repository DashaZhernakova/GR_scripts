#!/bin/bash

echo
echo "################" 
echo "Annotate GWAS results by adding variant-assoicated traits from GWAS catalog, ClinVar and HGMD"
echo "################"
echo

if [ $# -eq 0 ]
  then
    echo "No arguments supplied"
    echo "Usage: annotate_from_DBs.sh <input_file> <output_file>"
    echo "Exiting"
    exit 1
fi

inFile=$1
outFile=$2

echo "Annotating $inFile"
echo "Output will be in ${outFile}.txt"
echo

python ~/scripts/GWAS_scripts/gwas2qtlFile.py $inFile > ${inFile}.fake_qtl.txt

java -Xmx90g \
-jar /home/dzhernakova/tools/QtlAnnotator-1.0-SNAPSHOT-jar-with-dependencies.jar \
-a /home/dzhernakova/resources/variant_annotation/all_db_merged_31-07-18.codes.uniq.srt.hg19.txt \
-G VCF_FOLDER -g /home/dzhernakova/resources/b37/variation/1000G_biallelic_SNPs \
-q ${inFile}.fake_qtl.txt \
-o ${outFile}

rm ${inFile}.fake_qtl.txt
