#!/bin/bash

#
# Runs standard QC steps on genotype data: sex check, IBD, PCA, missing genotype rate, mendel error rate.
# Input is a VCF file; output is written to a new folder qc/ in the same folder as the VCF file
#

f=$1 # VCF file with genotypes to QC
folder=${f%"/"*}
mkdir ${folder}/qc/
fname=${f##*/}
out_base=${folder}/qc/${fname%.vcf.gz}
echo "Running genotype QC for $f"
echo "Output prefix: $out_base"

# Add family ids, split the parX regions
/mnt/sochi/dzhernakova/tools/plink_1.9/plink --vcf $f --split-x b38 --make-bed --out ${out_base}.splitx --update-ids /mnt/sochi/dzhernakova/GR/GR_freeze2019/sample_info/update_fid_iid.txt
# Add mother and father ids
plink --bfile ${out_base}.splitx --update-parents /mnt/sochi/dzhernakova/GR/GR_freeze2019/sample_info/update_parents.txt --update-sex /mnt/sochi/dzhernakova/GR/GR_freeze2019/sample_info/update_sex.txt --make-bed --out ${out_base}.splitx.pheno

# Check sex
/mnt/sochi/dzhernakova/tools/plink_1.9/plink --bfile ${out_base}.splitx.pheno --check-sex --out ${out_base}.splitx.pheno

# IBD
/mnt/sochi/dzhernakova/tools/plink_1.9/plink --bfile ${out_base}.splitx.pheno --maf 0.05 --geno 0.05 --hwe 1e-4 --make-bed -out ${out_base}.splitx.pheno.filtered
/mnt/pthome/dzhernakova/scripts/prune.sh ${out_base}.splitx.pheno.filtered
/mnt/sochi/dzhernakova/tools/plink_1.9/plink --bfile ${out_base}.splitx.pheno.filtered.pruned_r2_0.2 --genome --out ${out_base}.splitx.pheno.filtered.pruned

# PCA
/mnt/sochi/dzhernakova/tools/plink_1.9/plink --bfile ${out_base}.splitx.pheno --remove /mnt/sochi/dzhernakova/GR/GR_freeze2019/sample_info/children.txt --make-bed --out ${out_base}.no_children
plink --bfile ${out_base}.no_children --maf 0.05 --geno 0.05 --hwe 1e-4 --out ${out_base}.no_children.filtered --make-bed
/mnt/pthome/dzhernakova/scripts/prune.sh ${out_base}.no_children.filtered
plink --bfile ${out_base}.no_children.filtered.pruned_r2_0.2 --pca --out ${out_base}.no_children.filtered.pruned_r2_0.2

# missingness rate
plink --bfile ${out_base}.splitx.pheno --missing gz --out ${out_base}.splitx.pheno

# Mendel errors
plink --mendel summaries-only --bfile ${out_base}.splitx.pheno --out ${out_base}.splitx.pheno

# Draw plots
Rscript /mnt/sochi/dzhernakova/GR/GR_freeze2019/scripts/plot_all_QC.R ${out_base}