module load htslib
module load bcftools
module load gatk
module load java/1.8.0
module load plink

#chr=$1
chr=22

# Genotype calling
#/home/dzhernakova/GR/freeze2019/genotypes/recall_1000G_SNPs/call.sh \ # ../genotype_calling/call_at_positions.sh
#/home/dzhernakova/GR/freeze2019/genotypes/recall_1000G_SNPs/bam_list.txt \
#/home/dzhernakova/GR/freeze2019/genotypes/recall_1000G_SNPs//sampleids_gender.txt \
#chr${chr}.GR_bcftools_1kg_snps \
#/home/dzhernakova/GR/freeze2019/genotypes/recall_1000G_SNPs/ \
#1000000 5 \
#/home/dzhernakova/GR/freeze2019/genotypes/recall_1000G_SNPs/1kg_chr_pos/chr${chr}.chr_pos.tsv.gz \
#/home/dzhernakova/GR/freeze2019/genotypes/recall_1000G_SNPs/fasta_per_chr/chr${chr}.fa

in_dir=/home/dzhernakova/GR/freeze2019/genotypes/recall_1000G_SNPs/vcf/chr${chr}.GR_bcftools_1kg_snps/
cd $in_dir

# filtering 
java -jar ~/tools/picard-tools-2.5.0/picard.jar CreateSequenceDictionary R=../../fasta_per_chr/chr${chr}.fa O=../../fasta_per_chr/chr${chr}.dict
tabix -p vcf chr${chr}.GR_bcftools_1kg_snps.vcf.gz 

bcftools view -M2 -Oz -i 'QUAL>20' chr${chr}.GR_bcftools_1kg_snps.vcf.gz > chr${chr}.GR_bcftools_1kg_snps.QUAL20.vcf.gz
tabix -p vcf chr${chr}.GR_bcftools_1kg_snps.QUAL20.vcf.gz

/home/genomerussia/tools/gatk/gatk-4.1.0.0/gatk VariantFiltration \
-R ../../fasta_per_chr/chr${chr}.fa \
--variant chr${chr}.GR_bcftools_1kg_snps.QUAL20.vcf.gz \
-O chr${chr}.GR_bcftools_1kg_snps.QUAL20.DP5.vcf.gz \
--genotype-filter-expression "DP<5" \
--genotype-filter-name "DPFilter" \
--set-filtered-genotype-to-no-call true

bcftools annotate \
--rename-chrs ../../chr_renaming.txt \
-Oz chr${chr}.GR_bcftools_1kg_snps.QUAL20.DP5.vcf.gz \
> chr${chr}.GR_bcftools_1kg_snps.QUAL20.DP5.no_prefix.vcf.gz

tabix -p vcf chr${chr}.GR_bcftools_1kg_snps.QUAL20.DP5.no_prefix.vcf.gz



#
# Code below not tested!
#

 
# Convert vcf to plink 
plink \
 --vcf chr${chr}.GR_bcftools_1kg_snps.QUAL20.DP5.no_prefix.vcf.gz \
 --make-bed \
 --keep-allele-order \
 --update-ids ~/GR/freeze2019/sample_info/update_fid_iid.txt \
 --out tmp1.chr${chr}
 
#update phenotype information
plink \
 --bfile tmp1.chr${chr} \
 --make-bed \
 --update-parents ~/GR/freeze2019/sample_info/update_parents.txt \
 --update-sex ~/GR/freeze2019/sample_info/update_sex.txt \
 --exclude sample_to_exclude.txt \
 --keep-allele-order \
 --out tmp2.chr${chr}

# Filter on call rate, HWE and Mendel errors 
plink \
 --bfile tmp2.chr${chr} \
  --make-bed \
 --out chr${chr}.GR_bcftools_1kg_snps.QUAL20.DP5.no_prefix.geno0.05.hwe1e-4.me0.05 \
 --geno 0.05 --hwe 1e-4 --me 0.05 0.05 \
 --keep-allele-order
 
 rm tmp*.chr${chr}*
 
 # convert back to vcf
 plink \
 --bfile chr${chr}.GR_bcftools_1kg_snps.QUAL20.DP5.no_prefix.geno0.05.hwe1e-4.me0.05 \
 --recode vcf-iid bgz \
 --keep-allele-order \
 --out chr${chr}.GR_bcftools_1kg_snps.QUAL20.DP5.no_prefix.geno0.05.hwe1e-4.me0.05
 
 tabix -p vcf chr${chr}.GR_bcftools_1kg_snps.QUAL20.DP5.no_prefix.geno0.05.hwe1e-4.me0.05.vcf.gz
 

# Align GR genotypes to 1000G
java -Xmx50g \
-jar ~/tools/conform-gt.24May16.cee.jar \
ref=/home/dzhernakova/resources/1000G_b38/ALL.chr${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.SNPs.vcf.gz \
gt=chr${chr}.GR_bcftools_1kg_snps.QUAL20.DP5.no_prefix.geno0.05.hwe1e-4.me0.05.vcf.gz \
chrom=22 \
out=chr${chr}.GR_bcftools_1kg_snps.QUAL20.DP5.no_prefix.geno0.05.hwe1e-4.me0.05.conform1kg \
match=POS


# Phase
java -Xmx50g -jar /mnt/pthome/dzhernakova/tools/beagle.28Sep18.793.jar \
   ref=/home/dzhernakova/resources/1000G_b38/ALL.chr${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.SNPs.vcf.gz \
   gt=chr${chr}.GR_bcftools_1kg_snps.QUAL20.DP5.no_prefix.geno0.05.hwe1e-4.me0.05.conform1kg.vcf.gz \
   out=chr${chr}.GR_bcftools_1kg_snps.QUAL20.DP5.no_prefix.geno0.05.hwe1e-4.me0.05.beagle_phased \
   map=/mnt/pthome/dzhernakova/resources/recombination_maps/plink.chr${chr}.GRCh38.map \
   impute=false \
   nthreads=5
