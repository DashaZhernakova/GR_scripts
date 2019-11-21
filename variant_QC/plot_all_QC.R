args <- commandArgs(trailingOnly = TRUE)

#
# Plots QC metrics obtained by run_genotype_QC.sh script
#

file_base <- args[1] # file prefix of all QC results
region_file <- args[2] # file containg sample id in the first column and region in the second coumn
#file_base <- "/Users/dashazhernakova/Documents/Doby/GenomeRussia/freeze2019/genotype_QC/final2/GR_freeze2019.all_chr.QUAL40_GQ20_DP10_SP20.rsids.splitx.pheno"
#region_file <- "/Users/dashazhernakova/Documents/Doby/GenomeRussia/freeze2019/genotype_QC/final2/location-ethnicity.txt"
pdf(paste0(file_base, ".QC.pdf"))

#
# IBD
#
d = read.table(paste0(file_base, ".filtered.pruned.genome"), header = T)


par(pch=16, xpd=NA)
with(d,plot(Z0,Z1, xlim=c(0,1), ylim=c(0,1), type="n"))
with(subset(d,RT=="UN") , points(Z0,Z1, col = "grey", cex = 0.8))
with(subset(d,RT=="FS") , points(Z0,Z1,col=3))
with(subset(d,RT=="OT") , points(Z0,Z1,col=4, cex = 1.3))
with(subset(d,RT=="PO") , points(Z0,Z1,col=2, cex = 0.8))

legend(1,1, xjust=1, yjust=1, legend=c("full sibling", "other related", "parent-offspring", "unrelated"), pch=16, cex = 0.8, col=c(3,4,2,"grey"))

wrong1 <- d[((d$Z1< 0.7) & (d$RT == "PO")) | ((d$Z1 > 0.5) & (d$RT != "PO")),]
text(wrong1$Z0, wrong1$Z1, labels = wrong1$FID1, cex= 0.3, pos = 4)
print(wrong1)
#
# PCA
#

eigenvec <- read.table(paste0(file_base,".no_children.filtered.pruned_r2_0.2.eigenvec"), header = F)
colnames(eigenvec) = c(c("FID","IID"), paste0("PC", seq(5)))
eigenvec <- eigenvec[,seq(1,7)]
region_file <- "location-ethnicity.txt"
pop <- read.table(region_file, sep = "\t", header = T, as.is = T)

#pop_subset <- pop[match(mds_res$IID, pop$Sample_ID, nomatch = 0),]
pop_subset <- pop[match(eigenvec$IID, pop$Sample_ID, nomatch = 0),]
merged_tab <- merge(x=eigenvec, y=pop_subset, by.x = "IID", by.y = "Sample_ID")
merged_tab$Ethnicity <- as.factor(merged_tab$Ethnicity)


ggplot(merged_tab, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = Ethnicity, shape = Ethnicity), size=3)  + 
  xlab("PC1") + ylab("PC2")+ 
  theme(panel.background = element_rect(fill = "transparent")) + 
  scale_shape_manual( values=seq(0,17))




#merged_tab$Location <- as.factor(merged_tab$Ethnicity)
#pairs(merged_tab[,seq(3,7)], col = as.factor(merged_tab$Location), bg = rev(merged_tab$Location), pch = 21, upper.panel = NULL, cex = 0.8, labels = colnames(merged_tab)[seq(3,7)])
#legend("topright", bty = "n", legend = levels(merged_tab$Location), pch = 21, col=1:length(levels(merged_tab$Location)), xpd = T, cex = 0.8, y.intersp = 0.5)
#legend("topright", bty = "n", legend = levels(merged_tab$Location), pch = 21, col = rev(colrs), pt.bg=rev(colrs + 1), xpd = T, cex = 0.8, y.intersp = 0.5)
#legend("bottomright", bty = "n", legend = levels(merged_tab$Location), pch = 21, col=colrs, xpd = T, cex = 0.8, y.intersp = 0.5, bg=NA)

#
# Sex check
#

scheck <- read.table(paste0(file_base, ".splitx.pheno.sexcheck"), header = T)
plot(scheck$PEDSEX, scheck$SNPSEX, pch = 16, col = scheck$STATUS, xlab = "pheno sex", ylab = "geno sex", main = "Sex check")
legend("center", legend = c("Ok", "Problem"), pch = 16, col=c(1,2), xpd = T, cex = 0.8)

#
# Missing rate
#
miss <- read.tableread.table(paste0(file_base,".filtered.imiss"), header = T)
barplot(miss$F_MISS, ylim = c(0, max(0.1,max(miss$F_MISS))), names = miss$IID, ylab = "Missingness rate", main = "Fraction of missing genotypes after filtering", las=2, cex.names = 0.3, col = "steelblue1")
print(miss[miss$F_MISS > 0.1,])
#
# Mendel errors
#
mendel <- read.table(read.table(paste0(file_base,".filtered.fmendel"), header = T)
barplot(mendel$N, names = mendel$FID, ylab = "# mendel errors", las =3, cex.names = 0.4, col = "steelblue1", main = "# mendel errors in filtered genotypes")
print(mendel[mendel$N > 5000,])


dev.off()
