chr=$1

# bcftools view \
# -m2 -M2 \
# -Oz -o /mnt/sochi/dzhernakova/GR/GR_freeze2019/genotypes/per_chr/GR_freeze2019.filtered.biallelic.chr${chr}.vcf.gz  \
# /mnt/sochi/dzhernakova/GR/GR_freeze2019/genotypes/GR_freeze2019.filtered.vcf.gz  \
# "${chr}"
 
# tabix -p vcf /mnt/sochi/dzhernakova/GR/GR_freeze2019/genotypes/per_chr/GR_freeze2019.filtered.biallelic.chr${chr}.vcf.gz
 
# plink \
# --vcf /mnt/sochi/dzhernakova/GR/GR_freeze2019/genotypes/per_chr/GR_freeze2019.filtered.biallelic.chr${chr}.vcf.gz \
# --make-bed \
# --keep-allele-order \
# --update-ids /mnt/sochi/dzhernakova/GR/GR_freeze2019/sample_info/update_fid_iid.txt \
# --out tmp1.chr${chr}
 
# plink \
# --bfile tmp1.chr${chr} \
# --make-bed \
# --update-parents /mnt/sochi/dzhernakova/GR/GR_freeze2019/sample_info/update_parents.txt \
# --update-sex /mnt/sochi/dzhernakova/GR/GR_freeze2019/sample_info/update_sex.txt \
# --exclude sample_to_exclude.txt \
# --keep-allele-order \
# --out tmp2.chr${chr}
 
# plink \
# --bfile tmp2.chr${chr} \
#  --make-bed \
# --out GR_freeze2019.filtered.chr${chr}.maf0.01.geno0.05.hwe1e-4.me0.05 \
# --maf 0.01 --geno 0.05 --hwe 1e-4 --me 0.05 0.05 \
# --keep-allele-order
 
# rm tmp*.chr${chr}*
 
 plink \
 --bfile GR_freeze2019.filtered.chr${chr}.maf0.01.geno0.05.hwe1e-4.me0.05 \
 --recode vcf-iid bgz \
 --keep-allele-order \
 --out GR_freeze2019.filtered.chr${chr}.maf0.01.geno0.05.hwe1e-4.me0.05
 
 tabix -p vcf GR_freeze2019.filtered.chr${chr}.maf0.01.geno0.05.hwe1e-4.me0.05.vcf.gz
 
 java -Xmx50g -jar /mnt/pthome/dzhernakova/tools/beagle.28Sep18.793.jar \
   gt=GR_freeze2019.filtered.chr${chr}.maf0.01.geno0.05.hwe1e-4.me0.05.vcf.gz \
   out=GR_freeze2019.filtered.chr${chr}.maf0.01.geno0.05.hwe1e-4.me0.05.beagle_phased \
   map=/mnt/pthome/dzhernakova/resources/recombination_maps/plink.chr${chr}.GRCh38.map \
   impute=false \
   nthreads=5

tabix -p vcf GR_freeze2019.filtered.chr${chr}.maf0.01.geno0.05.hwe1e-4.me0.05.beagle_phased.vcf.gz

bcftools merge \
-Oz -o 1000G+GR.filtered.phased.chr${chr}.vcf.gz \
/mnt/sochi/dzhernakova/resources/1000G_b38/ALL.chr${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz \
GR_freeze2019.filtered.chr${chr}.maf0.01.geno0.05.hwe1e-4.me0.05.beagle_phased.vcf.gz \
> merge.chr${chr}.out


echo "Finished"

