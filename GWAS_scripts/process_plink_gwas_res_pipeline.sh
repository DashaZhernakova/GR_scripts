#!/bin/bash

# read file prefix (see below) of GWAS results from plink as the script 1st argument
# theresults will be written to $output_folder (2nd argument)

file_prefix=$1
output_folder=$2

if [ $# -ne 2 ]
  then
    echo "No arguments supplied"
    echo "Usage: process_plink_gwas_res_pipeline.sh <input_file_prefix> <output_folder>"
    echo "Exiting"
    exit 1
fi
mkdir $output_folder

#1. Add rs ids.
#/home/dzhernakova/scripts/GWAS_scripts/add_rsids.sh <file_prefix> <output_folder>
#
#where <file_prefix> is the path to the file + file name prefix. 
#The whole gwas files should conform to the following pattern:  ${file_prefix}.${test_name}.assoc.logistic.rsid.txt
#For example for the file 
#/home/dzhernakova/botswana/GWAS_res/wgs.array.imputed_covarWA/wgs.array.imputed_covarWA.recessive.assoc.logistic.rsid.txt 
#file_prefix=/home/dzhernakova/botswana/GWAS_res/wgs.array.imputed_covarWA/wgs.array.imputed_covarWA
#
#<output_folder> - folder to write the resuts
#

echo "Adding rs ids"
sh /home/dzhernakova/scripts/GWAS_scripts/add_rsids.sh $file_prefix $output_folder
returnCode=$?
if [ ! $returnCode -eq 0 ]; then
	echo "Adding rs ids failed, exiting"
	exit 1
fi

file_prefix_local=${output_folder}/${file_prefix##*/}
echo "file prefix for sorting the files: $file_prefix_local"

#2. Sort, cut and annotate GWAS results with gene names
#/home/dzhernakova/scripts/GWAS_scripts/sort_annotate_gwas_res.sh <file_prefix>
#

sh /home/dzhernakova/scripts/GWAS_scripts/sort_annotate_gwas_res.sh $file_prefix_local 
returnCode=$?
if [ ! $returnCode -eq 0 ]; then
	echo "Sorting failed, exiting"
	exit 1
fi

#3. Annotate with DBs
#/home/dzhernakova//scripts/GWAS_scripts/annotate_from_DBs.sh <input_file> <output_file>
#
#where <input_file> - plink GWAS result (possibly a subset of significant results annotated with genes or a file with all tests merged) from step 2
#<output_file> - file name of the annotated result


sh /home/dzhernakova//scripts/GWAS_scripts/annotate_from_DBs.sh  ${file_prefix_local}.all_tests.sorted.txt ${file_prefix_local}.all_tests.sorted.annot_pheno 
returnCode=$?
if [ ! $returnCode -eq 0 ]; then
	echo "Annotation failed, exiting"
	exit 1
fi
