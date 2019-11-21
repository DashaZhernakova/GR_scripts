library(ggplot2)
setwd("/Users/dashazhernakova/Documents/Doby/GenomeRussia/freeze2019/population_genetics/PCA/")
eigenvec <- read.table("all_chr.GR+egdp+sgdp.eurasian.no_children.filtered.pruned_r2_0.2.eigenvec", header = F)
colnames(eigenvec) = c(c("FID","IID"), paste0("PC", seq(5)))
eigenvec <- eigenvec[,seq(1,7)]

pop <- read.table("../GR+papers_samples.txt", sep = "\t", header = T, as.is = T)
pop_subset <- pop[match(eigenvec$IID, pop$IID, nomatch = 0),]
merged_tab <- merge(x=eigenvec, y=pop_subset, by.x = "IID", by.y = "IID")
merged_tab$population <- as.factor(merged_tab$population)
merged_tab$continent_region <- as.factor(merged_tab$continent_region)
merged_tab$source <- as.factor(merged_tab$source)


pdf("PC1-2_GR+papers_regions.pdf", width= 15, height = 10, useDingbats=FALSE)
ggplot(merged_tab, aes(x = PC1, y = PC2)) + geom_point(aes(color = continent_region), size=2)  + xlab("PC1") + ylab("PC2")
dev.off()

pdf("PC1-2_GR+papers_subset4_Europe.pdf", width= 15, height = 10, useDingbats=FALSE)
ggplot(merged_tab[merged_tab$continent_region=="Europe",], aes(x = PC1, y = PC2)) + geom_point(aes(color = population, shape = source), size=3)  + xlab("PC1") + ylab("PC2")+ theme(panel.background = element_rect(fill = "transparent")) 
dev.off()


pdf("PC1-2_GR+papers_Siberia.pdf", width= 15, height = 10, useDingbats=FALSE)
ggplot(merged_tab[merged_tab$continent_region=="Siberia",], aes(x = PC1, y = PC2)) + geom_point(aes(color = population, shape = source), size=3)  + xlab("PC1") + ylab("PC2")+ theme(panel.background = element_rect(fill = "transparent"))
dev.off()

pdf("PC1-2_GR.pdf", width= 15, height = 10, useDingbats=FALSE)
ggplot(merged_tab[merged_tab$source=="GR",], aes(x = PC1, y = PC2)) +
 geom_point(aes(color = population, shape = population), size=3)  +
 xlab("PC1") + ylab("PC2") +
 theme(panel.background = element_rect(fill = "transparent")) +
 scale_shape_manual(values=seq(0,43))
  
dev.off()


#
# batch effect test
#
pdf("batch_effect.pdf", width= 15, height = 10, useDingbats=FALSE)
pops <- c("Armenians", "Evens","Altaians",  "Adygei", "Khantys", "Yakuts", "Kyrgyz", "Uygur", "Komis", "Bashkirs", "Bengali", "Druze", "Brahmin")
ggplot(merged_tab[merged_tab$population %in% pops,], aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = population, shape = source), size=2)  + 
  xlab("PC1") + ylab("PC2")
dev.off()


#
# Europe
#
merged_tab2 <- merged_tab[merged_tab$PC1 > 0.02 & merged_tab$PC2 < 0,]
pdf("zoom_Europe_labels.pdf", width = 15, height = 15)
ggplot(merged_tab2, aes(x = PC1, y = PC2, label = population)) + 
 geom_point(aes(color = population), size=2) + geom_text() +
 xlab("PC1") + ylab("PC2") +
 scale_shape_manual(values=c(seq(0,25), seq(0,25), seq(0,25)))
dev.off()
#
# Look up
#
merged_tab[(merged_tab$PC2 > 0.038) & (merged_tab$continent_region == "Caucasus"), ]
query_pop = "Estonians"
ggplot(merged_tab[merged_tab$continent_region == "Europe",], aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = continent_region), size=1)  + 
  xlab("PC1") + ylab("PC2") + 
  geom_point(data = merged_tab[merged_tab$population == query_pop,], aes(x = PC1, y = PC2), colour="red", size=3)


#
# Plot differences between regions
#

pop <- read.table("../GR+papers_samples_3.txt", sep = "\t", header = T, as.is = T)
pop_subset <- pop[match(eigenvec$IID, pop$IID, nomatch = 0),]
merged_tab <- merge(x=eigenvec, y=pop_subset, by.x = "IID", by.y = "IID")
merged_tab$population <- as.factor(merged_tab$population)
merged_tab$continent_region <- as.factor(merged_tab$continent_region)
merged_tab$source <- as.factor(merged_tab$source)

#remove 2 ouliers from Komi:
merged_tab <- merged_tab[!merged_tab$IID %in% c("GS000035018-ASM", "GS000014328-ASM"),]
#rotate the PC space
merged_tab$PC1 <- -1*merged_tab$PC1
merged_tab$PC2 <- -1*merged_tab$PC2

europe <- merged_tab[merged_tab$PC1 > -0.04 & merged_tab$PC1 < 0 & merged_tab$PC2 > -0.02 & merged_tab$PC2 < 0.04,]


plots = list()
#par(mfrow = c(5,3))
for (sample_pop in unique(europe[europe$source == "GR","population"])){
  print (sample_pop)
  
  p <- ggplot(europe, aes(x = PC1, y = PC2)) + 
    geom_point(data = europe, color = "grey", size=1) +
    geom_point(data = europe[europe$population == sample_pop & europe$source == "GR",], aes(color = population_village), size=2) +
    xlab("PC1") + ylab("PC2") +
    scale_shape_manual(values=seq(0,43)) +
    ggtitle(sample_pop)
  plots <- c(plots, list(p))
  
}
pdf("PC1-2_GR+papers_europe_village.pdf", useDingbats=FALSE)
multiplot(plots)
dev.off()


siberia <-  merged_tab[merged_tab$PC1 > 0  & merged_tab$PC2 > 0,]


plots = list()
#par(mfrow = c(5,3))
for (sample_pop in unique(siberia[siberia$source == "GR","population"])){
  print (sample_pop)
  
  p <- ggplot(siberia, aes(x = PC1, y = PC2)) + 
    geom_point(data = siberia, color = "grey", size=1) +
    geom_point(data = siberia[siberia$population == sample_pop & siberia$source == "GR",], aes(color = population_village), size=2) +
    xlab("PC1") + ylab("PC2") +
    scale_shape_manual(values=seq(0,43)) +
    ggtitle(sample_pop)
  plots <- c(plots, list(p))
  
}
pdf("PC1-2_GR+papers_siberia_village.pdf", useDingbats=FALSE)
multiplot(plots)
dev.off()