#!/bin/bash

set -e
# List of BAM files
bam_list=$1

# Ped sample file
ped_file=$2

# Name prefix for the resulting VCF
vcf_prefix=$3

# Results path
results_path=$4

# Region size for mpileup
mpileup_region_size=$5

# bcftools thread number
bcftools_threads=$6

# regions to call variants in
loci_to_call=$7
# Tools and parameters

tools_path="/mnt/gr1/tools/"
samtools="$tools_path/samtools-1.9/samtools"
bcftools_dir="$tools_path/bcftools-1.9/"
bcftools="$tools_path/bcftools-1.9/bcftools"
mpileup_Q="20"
mpileup_q="20"
samtools_field_list="AD,INFO/AD,DP,SP"
bcftools_field_list="GQ,GP"
min_QUAL="40"
min_GQ="20"
min_DP="10"
max_SP="20"

# Data

ref_fa="/mnt/gr1/fasta/Homo_sapiens.GRCh38.dna.toplevel.main-chr.chr.fa"
ref_fai="/mnt/gr1/fasta/Homo_sapiens.GRCh38.dna.toplevel.main-chr.chr.fa.fai"
#ref_fa="/mnt/gr1/genotyping/hg38.chr10.10mb.fa"
#ref_fai="/mnt/gr1/genotyping/hg38.chr10.10mb.fa.fai"
#loci_to_call="/mnt/gr1/resources/gatk_intervals.bed"

# Generate a random ID for this run
random_id=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 15 | head -1`

log="$results_path/log"
if [ ! -d "$log/log-$random_id" ]; then
    mkdir -p $log/log-$random_id
fi

tmp="$results_path/tmp"
if [ ! -d "$tmp/tmp-$random_id" ]; then
    mkdir -p $tmp/tmp-$random_id
fi

vcf="$results_path/vcf"
if [ ! -d "$vcf/$vcf_prefix" ]; then
    mkdir -p $vcf/$vcf_prefix
fi

# Print input info

date
echo
echo "Pipeline starts with the following parameters:"
echo "----------------------------------------------"
echo "BAM list:                $bam_list"
echo "VCF file prefix:         $vcf_prefix"
echo "----------------------------------------------"
echo "Region size for mpileup: $mpileup_region_size"
echo "bcftools thread number:  $bcftools_threads"
echo "ped file: $ped_file"
echo "Minimum base quality:    $mpileup_Q"
echo "Samtools mpileup fields: $samtools_field_list"
echo "Bcftools fields:         $bcftools_field_list"
echo "min QUAL:                $min_QUAL"
echo "min GQ:                  $min_GQ"
echo "min DP:                  $min_DP"
echo "max SP:                  $max_SP"
echo "Reference FASTA:         $ref_fa"
echo "Reference FAI:           $ref_fai"
echo "log directory:           $log/log-$random_id"
echo "tmp directory:           $tmp/tmp-$random_id"
echo "vcf directory:           $vcf/$vcf_prefix"
echo


# Call SNPs
echo "-> Call SNPs..."


time ${tools_path}/vcfutils.pl splitchr -l $mpileup_region_size $ref_fai | \
    xargs -I {} -P $bcftools_threads sh -c \
    "$bcftools mpileup \
	-b $bam_list \
        -a $samtools_field_list \
        -Q $mpileup_Q \
        -q $mpileup_q \
	-f $ref_fa \
        -Ou \
        -r '{}' \
        -T ${loci_to_call} \
        --skip-indels \
        2> $log/log-$random_id/\"err.samtools.mpileup.{}.log\" | \
    $bcftools call \
        -Ou \
        -v \
        -m \
        -f $bcftools_field_list \
        --ploidy GRCh38 \
        -V indels \
        -S ${ped_file} \
        > $tmp/tmp-$random_id/\"tmp.{}.bcf\" \
        2> $log/log-$random_id/\"err.bcftools.call.{}.log\""

returnCode=$?
echo "bcftools return code: $returnCode"
if [ $returnCode -ne 0 ]
then
	#Return non zero return code
	exit 1
fi

echo "Done."
echo

# Concatenate VCF files from all partial tmp.*.bcf files, rename samples
echo "-> Concatenate VCF files from all partial tmp.*.bcf files and rename samples in the concatenated VCF file..."


cat $bam_list | \
while read line; do
  fname=${line##*/}
  sample_id=${fname%.rmdup.bam}
  echo "$line $sample_id" >> $tmp/tmp-$random_id/sample_renaming.txt
done

time $bcftools concat \
  -O u \
  --threads $bcftools_threads \
  `ls $tmp/tmp-$random_id/tmp.*.bcf | sort -V` | \
$bcftools reheader \
    -s $tmp/tmp-$random_id/sample_renaming.txt \
    - | \
    $bcftools view \
    -O z \
    -o $vcf/$vcf_prefix/$vcf_prefix.vcf.gz \
    - \
2> $log/log-$random_id/err.concatenate_and_rename_samples_vcfs.log

echo "Done."
echo

# Remove all temporary bcf files
echo "-> Clean $tmp/tmp-$random_id/tmp.*.bcf..."
#rm $tmp/tmp-$random_id/tmp.*.bcf
#rm $tmp/tmp-$random_id/sample_renaming.txt
echo "Done."
echo

# Filter VCF file
echo "-> Filter VCF file..."
time $bcftools view -i "%QUAL>$min_QUAL && FORMAT/GQ>$min_GQ && FORMAT/DP>$min_DP && FORMAT/SP<$max_SP" \
                    -O z \
                    -o $vcf/$vcf_prefix/$vcf_prefix.sorted.minQUAL"$min_QUAL".minGQ"$min_GQ".minDP"$min_DP".maxSP"$max_SP".vcf.gz \
                    $vcf/$vcf_prefix/$vcf_prefix.vcf.gz \
                    2> $log/log-$random_id/err.bcftools.view.log
echo "Done."
echo

echo "Finished."
date
 