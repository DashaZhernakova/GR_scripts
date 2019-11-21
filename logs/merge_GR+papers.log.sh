#/bin/bash

#
# Merge GR genotypes with Mallick and Pagani, apply basic filtering
#

for chr in `seq 1 22`
do
 
 # merge GR + Pagani et al
 /mnt/pthome/dzhernakova/scripts/merge_excl_both.sh \
 ../per_chr/GR_freeze2019.filtered.pheno.chr${chr}  \
 /mnt/sochi/dzhernakova/GR/from_papers/pagani/egdp_remapped/b38/per_chr/multi-sample.b38.chr${chr} \
 chr${chr}.GR+egdp \
 > merge.chr${chr}.GR+egdp.log 2>&1
 
 # merge GR and Pagani et al + Mallick et al.
 /mnt/pthome/dzhernakova/scripts/merge_excl_both.sh \
 chr${chr}.GR+egdp \
 /mnt/sochi/dzhernakova/GR/from_papers/simons/sgdp_variants/b38/cteam_extended.v4.PS3_phase.public.b38.sorted.rmdup.chr${chr} \
 chr${chr}.GR+egdp+sgdp \
 > merge.chr${chr}.GR+egdp+sgdp.log 2>&1

 echo "chr${chr}.GR+egdp+sgdp" >> flist.txt

done

plink --merge-list flist.txt --make-bed --out all_chr.GR+egdp+sgdp

# exclude low mappability regions (Umap mask)
sh /mnt/pthome/dzhernakova/scripts/extract_bed_from_plink.sh all_chr.GR+egdp+sgdp /mnt/gr1/resources/k100.umap.bed.gz all_chr.GR+egdp+sgdp.mapp

#exclude non-variant sites
awk '{if (($5==0) || ($6==0)) print $2}' all_chr.GR+egdp+sgdp.mapp.bim > non_variant_sites.txt
plink --bfile all_chr.GR+egdp+sgdp.mapp --exclude non_variant_sites.txt --make-bed --out all_chr.GR+egdp+sgdp.mapp.vars