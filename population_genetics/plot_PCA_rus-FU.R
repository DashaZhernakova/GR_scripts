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


#remove 2 ouliers from Komi:
merged_tab <- merged_tab[!merged_tab$IID %in% c("GS000035018-ASM", "GS000014328-ASM"),]
#rotate the PC space
merged_tab$PC1 <- -1*merged_tab$PC1
merged_tab$PC2 <- -1*merged_tab$PC2

pdf("PC1-2_GR+papers_eurasian_region_rotated.pdf", width= 15, height = 10, useDingbats=FALSE)
ggplot(merged_tab, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = continent_region), size=2) +
  xlab("PC1") + ylab("PC2") +
  scale_shape_manual(values=seq(0,25))

dev.off()


pdf("PC1-2_GR+papers_finno-ugric-russian.pdf", width= 15, height = 10, useDingbats=FALSE)
ggplot(merged_tab, aes(x = PC1, y = PC2)) + 
  geom_point(color = "grey", size=1) + 
  geom_point(data = merged_tab[merged_tab$grouping != "0",], aes(color = population, shape = population), size=2) +
  xlab("PC1") + ylab("PC2") +
  scale_shape_manual(values=seq(0,25))

dev.off()

merged_tab$grouping <- as.factor(merged_tab$grouping)
cbp1 <- c("#E69F00", "#56B4E9", "#009E73",
          "#FF0000", "#0072B2", "#D55E00", "#CC79A7")
pdf("PC1-2_GR+papers_finno-ugric-russian2.pdf", width= 15, height = 10, useDingbats=FALSE)
ggplot(merged_tab, aes(x = PC1, y = PC2)) + 
  geom_point(color = "grey", size=1) + 
  geom_point(data = merged_tab[merged_tab$grouping2 != "0",], aes(color = grouping2, shape = grouping2), size=3) +
  xlab("PC1") + ylab("PC2") +
  scale_shape_manual(values=c(16,16,16,17,16,16,16)) +
  scale_color_manual(values=cbp1) +
  theme(legend.text=element_text(size=12), legend.title = element_blank())
dev.off()


merged_tab[merged_tab$population %in% c("Khantys", "Mansis", "Seto", "Hungarians", "Saami", "Maris", "Udmurts", "Komis"),]$grouping <- "FU2"
merged_tab$grouping <- as.factor(merged_tab$grouping)

myPalette <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#9900FF",  
               "#FF0000", "#D55E00", "#CC79A7", "#CC0000", "#00CCFF",  "#FF9933", "#FF00CC", "#990000", "#663300")
# "#66FFFF", "#000999", "#009966", "#000066",
myPalette <- c("#E69F00", "#00CC99", "#009E73", "#66FF00", "#00FFCC", "#CC00CC", "#3333CC", "#9933FF", "#CC6633", 
               "#666633",
               "#CC3300", "#000066", "#FF3333",  "#FF9933")
myShapes <- c(0, 16, 17, 18, 16, seq(1, 4), 15, seq(5,7))
pdf("PC1-2_GR+papers_finno-ugric-russian_zoom2.pdf", width= 15, height = 10, useDingbats=FALSE)
ggplot(merged_tab, aes(x = PC1, y = PC2)) + 
  geom_point(data = merged_tab[!merged_tab$grouping %in% c("Rus", "FU"),], color = "grey", size=1) + 
  geom_point(data = merged_tab[merged_tab$grouping %in% c("Rus", "FU"),], aes(color = population, shape = population), size=5) +
  xlab("PC1") + ylab("PC2") + xlim(-0.04, -0.02) + ylim(0.004,0.02) +
  scale_shape_manual(values=myShapes) +
  scale_color_manual(values=myPalette)

dev.off()

pdf("PC1-2_GR+papers_finno-ugric-russian_zoom.pdf", width= 15, height = 10, useDingbats=FALSE)
ggplot(merged_tab, aes(x = PC1, y = PC2)) + 
  geom_point(data = merged_tab[!merged_tab$grouping %in% c("Rus", "FU"),], color = "grey", size=1) + 
  geom_point(data = merged_tab[merged_tab$grouping %in% c("Rus", "FU"),], aes(color = population, shape = population), size=3) +
  xlab("PC1") + ylab("PC2") + xlim(-0.04, 0) + ylim(0.004,0.04) +
  scale_shape_manual(values=seq(0,43)) +

dev.off()

eur_subs <- merged_tab[merged_tab$PC1 > -0.04 & merged_tab$PC1 < 0 & merged_tab$PC2 > 0 & merged_tab$PC2 < 0.04,]
pdf("PC1-2_GR+papers_finno-ugric-russian_zoom_all_pop.pdf", width= 15, height = 10, useDingbats=FALSE)
ggplot(eur_subs, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = population, shape = population), size=3) +
  xlab("PC1") + ylab("PC2") + xlim(-0.04, 0) + ylim(0.004,0.04) +
  scale_shape_manual(values=seq(0,43))

dev.off()

fu <- c("Estonians", "Finnish", "Ingrians", "Karelians", "Khantys", "Komis", "Mansis", "Maris", "Saami",  "Udmurts", "Vepsas", "Hungarians")
fu_subs <- c("Estonians", "Finnish", "Ingrians", "Karelians",  "Komis", "Maris", "Saami",  "Udmurts", "Vepsas")
sam <- c("Forest-Nenets", "Nganasans", "Selkups", "Tundra-Nenets")
rus <- c("Arkhangelsk", "Novgorod", "Pskov", "Rostov", "Rus_Udmurtia", "Vladimir", "Voronezh", "Yaroslavl")

