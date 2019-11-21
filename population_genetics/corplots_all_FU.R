library("reshape2")
library(corrplot)
setwd("/Users/dashazhernakova/Documents/Doby/GenomeRussia/freeze2019/population_genetics/Dstats/")
d <- read.table("fu_vs_fu.txt", header = T, sep = "\t", as.is = T, check.names = F)

#pop_order <- c("Hungarians", "Estonians", "Seto", "Finnish", "Vepsas", "Ingrians", "Karelians", "Saami", "Mansis", "Maris", "Udmurts", "Komis", "Khantys", "Tundra-Nenets", "Forest-Nenets", "Selkups", "Nganasans")
pop_order <- c("Hungarians", "Estonians", "Seto", "Finnish", "Vepsas", "Ingrians", "Karelians", "Saami", "Mordvins","Maris", "Udmurts", "Komis", "Mansis", "Khantys")
rus <- c("Novgorod", "Pskov", "Arkhangelsk", "Rus_Udmurtia", "Yaroslavl", "Vladimir", "Voronezh", "Rostov")

pdf("all_rus_vs_FU_Dstats.pdf", width = 10, height = 30, useDingbats=F)
par(mfrow=c(4,2))
for (pop in rus){
  print(pop)
 z_matr <- acast(d[d$P2 == pop,], P3~P4, value.var = "z")
 
 z_matr[is.na(z_matr)] = 0
 z_matr[which(z_matr>-3 & z_matr < 3)]=0
 z_matr[which(z_matr > 15)] = 15
 z_matr[which(z_matr < -15)] = -15
 
 z_matr <- z_matr[match(pop_order, row.names(z_matr), nomatch = 0), match(pop_order, colnames(z_matr), nomatch = 0)]
 corrplot(z_matr, tl.cex=1, tl.col=gray(0.3), 
         cl.lim=c(-15, 15), is.corr=F, rect.lwd=0.2,  tl.offset = 0.5, 
         cl.cex = 0.8, cl.ratio = 0.3, type = "upper", title = pop, mar=c(0,0,1,0))
}
dev.off()