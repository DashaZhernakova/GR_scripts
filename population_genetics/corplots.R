library(corrplot)
library(reshape2)
setwd("/Users/dashazhernakova/Documents/Doby/GenomeRussia/freeze2019/population_genetics/Dstats/")
col1<-colorRampPalette(c( "darkred", "#FF7F00", "white",  "#007FFF", "darkblue"))

fu_order <- c("Latvians", "Lithuanians","Estonians", "Finnish", "Vepsas", "Ingrians", "Karelians", "Saami", "Mordvins","Maris","Udmurts", "Komis",  "Mansis", "Khantys")
rus_order <- c("Yaroslavl", "Vladimir", "Voronezh", "Rostov", "Novgorod", "Pskov", "Arkhangelsk", "Rus_Udmurtia", "Udmurts")


d <- read.table("rostov_vs_fu.txt", header = T, sep = "\t", as.is = T, check.names = F)
pdf("corrplot_rostov_fu.pdf", height = 10, width = 15, useDingbats=F)
par(mfrow=c(1,2))
z_matr <- acast(d, P2~P4, value.var = "z")

z_matr <- z_matr[match(fu_order, row.names(z_matr), nomatch = 0), match(rus_order, colnames(z_matr), nomatch = 0)]

z_matr[which(z_matr > -3 & z_matr < 3)] = 0
z_matr[is.na(z_matr)] = 0
z_matr[which(z_matr > 15)] = 15
z_matr[which(z_matr < -15)] = -15
corrplot(z_matr, tl.cex=1, tl.col=gray(0.3), type='full', 
         cl.lim=c(-15, 15), is.corr=F, rect.lwd=0.2,  tl.offset = 0.5, 
         cl.cex = 0.8, cl.ratio = 0.3)


z_matr <- acast(d, P2~P4, value.var = "z")
z_matr <- z_matr[match(fu_order, row.names(z_matr), nomatch = 0), match(rus_order, colnames(z_matr), nomatch = 0)]
z_matr[is.na(z_matr)] = 0

corrplot(z_matr, tl.cex=1, tl.col=gray(0.3), type='full', 
         cl.lim=c(min(z_matr), max(z_matr)), is.corr=F, rect.lwd=0.2,  tl.offset = 0.5, 
         cl.cex = 0.8, cl.ratio = 0.3, method = "number")
mtext("D(Mbuti, X (row; Rostov, RUS (column) )", side = 3, line = -2, outer = TRUE, cex = 2)
dev.off()


d <- read.table("voronezh_vs_fu.txt", header = T, sep = "\t", as.is = T, check.names = F)
pdf("corrplot_voronezh_fu2.pdf", height = 10, width = 15, useDingbats=F)
par(mfrow=c(1,2))
z_matr <- acast(d, P2~P4, value.var = "z")

z_matr <- z_matr[match(fu_order, row.names(z_matr), nomatch = 0), match(rus_order, colnames(z_matr), nomatch = 0)]

z_matr[which(z_matr > -2.5 & z_matr < 2.5)] = 0
z_matr[is.na(z_matr)] = 0
z_matr[which(z_matr > 15)] = 15
z_matr[which(z_matr < -15)] = -15
corrplot(z_matr, tl.cex=1, tl.col=gray(0.3), type='full', 
         cl.lim=c(-15, 15), is.corr=F, rect.lwd=0.2,  tl.offset = 0.5, 
         cl.cex = 0.8, cl.ratio = 0.3)


z_matr <- acast(d, P2~P4, value.var = "z")
z_matr <- z_matr[match(fu_order, row.names(z_matr), nomatch = 0), match(rus_order, colnames(z_matr), nomatch = 0)]
z_matr[is.na(z_matr)] = 0
z_matr[which(z_matr > 15)] = 15
z_matr[which(z_matr < -15)] = -15
corrplot(z_matr, tl.cex=1, tl.col=gray(0.3), type='full', 
         cl.lim=c(-15,15), is.corr=F, rect.lwd=0.2,  tl.offset = 0.5, 
         cl.cex = 0.8, cl.ratio = 0.3, method = "number")
mtext("D(Mbuti, X (row); Voronezh, RUS (column) )", side = 3, line = -2, outer = TRUE, cex = 2)
dev.off()













pdf("corrplot_FU.pdf", height = 10, width = 15, useDingbats=F)
par(mfrow=c(1,2))
d <- as.matrix(read.table("dstats_fu.txt", header = T, sep = "\t", as.is = T, check.names = F, row.names = 1))
d[which(d>-3 & d < 3)]=0

corrplot(d, tl.cex=1, tl.col=gray(0.3), type='full', 
         cl.lim=c(-10, 10), is.corr=F, rect.lwd=0.2,  tl.offset = 0.5, 
         cl.cex = 0.8, cl.ratio = 0.3)

d <- as.matrix(read.table("dstats_fu.txt", header = T, sep = "\t", as.is = T, check.names = F, row.names = 1))
d[which(d>-2.5 & d < 2.5)]=0

corrplot(d, tl.cex=1, tl.col=gray(0.3), type='full', 
         cl.lim=c(-10, 10), is.corr=F, rect.lwd=0.2,  tl.offset = 0.5, 
         cl.cex = 0.8, cl.ratio = 0.3, method = "number")

dev.off()


pdf("corrplot_rus_arkhangelsk.pdf", height = 10, width = 15, useDingbats=F)
par(mfrow=c(1,2))
d <- as.matrix(read.table("dstats_rus_udm3.txt", header = T, sep = "\t", as.is = T, check.names = F, row.names = 1))
d[which(d>-3 & d < 3)]=0
d[which(d > 15)] = 15

corrplot(d, tl.cex=1, tl.col=gray(0.3), type='full', 
  cl.lim=c(-15, 15), is.corr=F, rect.lwd=0.2,  tl.offset = 0.5, 
  cl.cex = 0.8, cl.ratio = 0.3)

d <- as.matrix(read.table("dstats_rus_udm3.txt", header = T, sep = "\t", as.is = T, check.names = F, row.names = 1))
d[which(d>-2.5 & d < 2.5)]=0
d[which(d > 15)] = 15

corrplot(d, tl.cex=1, tl.col=gray(0.3), type='full', 
         cl.lim=c(-15, 15), is.corr=F, rect.lwd=0.2,  tl.offset = 0.5, 
         cl.cex = 0.8, cl.ratio = 0.3, method = "number")

dev.off()


d <- read.table("fu_west_east.txt", header = T, sep = "\t", as.is = T, check.names = F)
d <- d[d$P3 == "Voronezh",]
pdf("fu_west_east2.pdf", height = 10, width = 18, useDingbats=F)
par(mfrow=c(1,2))
z_matr <- acast(d, P2~P4, value.var = "z")

z_matr <- z_matr[match(fu_order, row.names(z_matr), nomatch = 0), match(fu_order, colnames(z_matr), nomatch = 0)]

z_matr[which(z_matr > -3 & z_matr < 3)] = 0
z_matr[is.na(z_matr)] = 0
z_matr[which(z_matr > 15)] = 15
z_matr[which(z_matr < -15)] = -15
corrplot(z_matr, tl.cex=1, tl.col=gray(0.3), type='full', 
         cl.lim=c(-15, 15), is.corr=F, rect.lwd=0.2,  tl.offset = 0.5, 
         cl.cex = 0.8, cl.ratio = 0.3, mar=c(0,0,1,0))

z_matr <- acast(d, P2~P4, value.var = "z")
z_matr <- z_matr[match(fu_order, row.names(z_matr), nomatch = 0), match(fu_order, colnames(z_matr), nomatch = 0)]
z_matr[is.na(z_matr)] = 0
z_matr[which(z_matr > 15)] = 15
z_matr[which(z_matr < -15)] = -15
corrplot(z_matr, tl.cex=1, tl.col=gray(0.3), type='full', 
         cl.lim=c(-15,15), is.corr=F, rect.lwd=0.2,  tl.offset = 0.5, 
         cl.cex = 0.8, cl.ratio = 0.3, mar=c(0,0,1,0), method = "number")
mtext("D(Mbuti, FU row; Voronezh, FU column )", side = 3, line = -2, outer = TRUE, cex = 2)
dev.off()
