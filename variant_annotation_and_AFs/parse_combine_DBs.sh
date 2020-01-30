#/bin/bash

#
# Downloads GWAS catalog, (HGMD) and ClinVar, reformats them in a unified way, merges annotations for the same SNPs
#

# GWAS catalog
wget https://www.ebi.ac.uk/gwas/api/search/downloads/alternative
cut -f12,13,22,35 alternative | grep -v " x " |  awk 'BEGIN {FS=OFS="\t"}; {if ((length($1) > 1) && (length($4) > 1)) print}' | tail -n+2 | sort -k1,1 -k2,2n > gwas_catalog_v1.0.2-associations_e96_r2019-05-03_cut.txt
python ../split_gwascatalog_multisnp_rows.py gwas_catalog_v1.0.2-associations_e96_r2019-05-03_cut.txt  > gwas_catalog_v1.0.2-associations_e96_r2019-05-03_cut_fmt.txt
python ../merge_gwascatalog_snps.py gwas_catalog_v1.0.2-associations_e96_r2019-05-03_cut_fmt.txt > gwas_catalog_v1.0.2-associations_e96_r2019-05-03_cut_fmt_merged.txt

# HGMD
#bcftools query -i 'TYPE="snp"' -f '%CHROM\t%POS\t%INFO/DB\t%INFO/CLASS\:%INFO/PHEN\n' /home/gtamazian/hgmd_pro_2018.2_hg38.vcf | awk 'BEGIN {FS=OFS="\t"}; {if ($3==".") $3 = $1 ":" $2; print}' > hgmd_cut_23-05-19.txt


# ClinVar
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
zgrep "GRCh38" variant_summary.txt.gz | grep "single nucleotide variant" | awk 'BEGIN {FS=OFS="\t"}; {print $19, $20, "rs" $10, $14 " (" $13 ")"}' | awk 'BEGIN {FS=OFS="\t"}; {if ($3=="rs-1") $3 = $1 ":" $2; print }' | grep -v "not specified" > clinvar_snps_cut_230519.txt
python ../merge_gwascatalog_snps.py  clinvar_snps_cut_230519.txt >  clinvar_snps_cut_230519_merged.txt 


awk 'BEGIN {FS=OFS="\t"}; {$4 = "GW:" $4; print}' gwas_catalog_v1.0.2-associations_e96_r2019-05-03_cut_fmt_merged.txt | sed 's/;/;GW:/g' > gwas_catalog_v1.0.2-associations_e96_r2019-05-03_cut_fmt_merged.GW_code.txt
awk 'BEGIN {FS=OFS="\t"}; {$4 = "CV:" $4; print }' clinvar_snps_cut_230519_merged.txt > clinvar_snps_cut_230519_merged.CV_code.txt
#awk 'BEGIN {FS=OFS="\t"}; {$4 = "HGMD:" $4; print }' hgmd_cut_23-05-19.txt > hgmd_cut_23-05-19.HGMD_code.txt

echo -e "chr\tpos\tid\ttrait" > all_db_merged_23-05-19.codes.txt
cat gwas_catalog_v1.0.2-associations_e96_r2019-05-03_cut_fmt_merged.GW_code.txt >> all_db_merged_23-05-19.codes.txt
cat clinvar_snps_cut_230519_merged.CV_code.txt >> all_db_merged_23-05-19.codes.txt
#cat hgmd_cut_23-05-19.HGMD_code.txt >> all_db_merged_23-05-19.codes.txt
python ../merge_gwascatalog_snps.py all_db_merged_23-05-19.codes.txt > all_db_merged_23-05-19.codes.uniq.txt

head -1 all_db_merged_23-05-19.codes.uniq.txt > all_db_merged_23-05-19.codes.uniq.srt.txt
tail -n+2 all_db_merged_23-05-19.codes.uniq.txt | sort -k1,1 -k2,2n >> all_db_merged_23-05-19.codes.uniq.srt.txt
bgzip all_db_merged_23-05-19.codes.uniq.srt.txt
tabix -f -b 2 -s 1 -e 2 -S1 all_db_merged_23-05-19.codes.uniq.srt.txt.gz 