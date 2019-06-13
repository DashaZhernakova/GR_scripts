library(gdsfmt)
#dat.f<-openfn.gds("/Users/malovs/Desktop/GWAS/RG/DATA/SPRG_05_02.gds")
dat.f<-openfn.gds("/mnt/genomerussia/dzhernakova/GR/ps+nov+yak+papers.filtered.CR0.05.MAF0.05.HWE1e-4.pruned_r2_0.2.gds")
gen<-read.gdsn(index.gdsn(dat.f,"genotype"))
ant<-read.gdsn(index.gdsn(dat.f,"sample.annot"))
closefn.gds(dat.f)
sors<-ant$Source
#
chisq_tst<-function(genotype,phtpe)
{
num.g<-dim(genotype)[2]	
qas<-NULL
pv<-NULL
for (i in 1:num.g) 
{
 		gtype<-genotype[,i]    
		sbst<-(gtype!=3)  
		gtype<-gtype[sbst]
		phtype<-phtpe[sbst]	
 		c.tab<-table(phtype,gtype)
 		if (min(dim(c.tab))>1){
			sr <- rowSums(c.tab)
	    	sc <- colSums(c.tab)
	    	E <- outer(sr, sc, "*")/sum(sr)    
			pv[i]<-chisq.test(phtype,gtype)$p.value 		
		} else pv[i]<-NA 
}
data.out<-list()
data.out$pv<-pv
} 
# ALL
pv<-chisq_tst(gen,sors)
#PV<-data.frame(cd=pv)
# SP
ftr<-sors!="GR"
pvsp<-chisq_tst(gen[ftr,],sors[ftr])
#
dat.f<-openfn.gds("/Users/malovs/Desktop/GWAS/RG/DATA/SPRG_05_02.gds",readonly=FALSE)
id<-read.gdsn(index.gdsn(dat.f,"sample.id"))
snp<-read.gdsn(index.gdsn(dat.f,"snp.id"))
chr<-read.gdsn(index.gdsn(dat.f,"snp.chromosome"))
pos<-read.gdsn(index.gdsn(dat.f,"snp.position"))
all<-read.gdsn(index.gdsn(dat.f,"snp.allele"))
gen<-read.gdsn(index.gdsn(dat.f,"genotype"))
ant<-read.gdsn(index.gdsn(dat.f,"sample.annot"))
add.gdsn(dat.f,"pv.source", val=pv,compress="ZIP",closezip=TRUE)
add.gdsn(dat.f,"pv.sp", val=pvsp,compress="ZIP",closezip=TRUE)
closefn.gds(dat.f)	


#
dat.f<-openfn.gds("/Users/malovs/Desktop/GWAS/RG/DATA/SPRG_05_02.gds",readonly=FALSE)
id<-read.gdsn(index.gdsn(dat.f,"sample.id"))
snp<-read.gdsn(index.gdsn(dat.f,"snp.id"))
chr<-read.gdsn(index.gdsn(dat.f,"snp.chromosome"))
pos<-read.gdsn(index.gdsn(dat.f,"snp.position"))
all<-read.gdsn(index.gdsn(dat.f,"snp.allele"))
gen<-read.gdsn(index.gdsn(dat.f,"genotype"))
ant<-read.gdsn(index.gdsn(dat.f,"sample.annot"))
pv<-read.gdsn(index.gdsn(dat.f,"pv.sp"))
closefn.gds(dat.f)	
#
pvb<-0.05/length(snp)
ftr<-pv>=pvb
snp<-snp[ftr]
chr<-chr[ftr]
pos<-pos[ftr]
all<-all[ftr]
gen<-gen[,ftr]
#
dat.f<-createfn.gds("/Users/malovs/Desktop/GWAS/RG/DATA/SPRG_05_02_pvcor.gds")
add.gdsn(dat.f,"sample.id", val=id,compress="ZIP",closezip=TRUE)
add.gdsn(dat.f,"snp.id", val=snp,compress="ZIP",closezip=TRUE)
add.gdsn(dat.f,"snp.chromosome", val=chr,compress="ZIP",closezip=TRUE)
add.gdsn(dat.f,"snp.position", val=pos,compress="ZIP",closezip=TRUE)
add.gdsn(dat.f,"snp.allele", val=all,compress="ZIP",closezip=TRUE)
add.gdsn(dat.f,"genotype", val=gen,compress="ZIP",closezip=TRUE,storage="bit2",replace=TRUE)
add.gdsn(dat.f,"sample.annot", val=ant,compress="ZIP",closezip=TRUE)
put.attr.gdsn(index.gdsn(dat.f,"genotype"),"sample.order")
closefn.gds(dat.f)	
# sort
dat.f<-openfn.gds("/Users/malovs/Desktop/GWAS/RG/DATA/SPRG_05_02_pvcor.gds",readonly=FALSE)
id<-read.gdsn(index.gdsn(dat.f,"sample.id"))
snp<-read.gdsn(index.gdsn(dat.f,"snp.id"))
chr<-read.gdsn(index.gdsn(dat.f,"snp.chromosome"))
pos<-read.gdsn(index.gdsn(dat.f,"snp.position"))
all<-read.gdsn(index.gdsn(dat.f,"snp.allele"))
gen<-read.gdsn(index.gdsn(dat.f,"genotype"))
ant<-read.gdsn(index.gdsn(dat.f,"sample.annot"))
closefn.gds(dat.f)	
# 
ant$ord<-0
ant$ord[ant$Population=="Yakut"]<-1
ant$ord[ant$Population=="Pskov"]<-2
ant$ord[ant$Population=="Novgorod"]<-3
ftr<-order(ant$ord)
ant<-ant[ftr,]
id<-id[ftr]
gen<-gen[ftr,]
dat.f<-createfn.gds("/Users/malovs/Desktop/GWAS/RG/DATA/SPRG_05_02_pvcor_ordered.gds")
add.gdsn(dat.f,"sample.id", val=id,compress="ZIP",closezip=TRUE)
add.gdsn(dat.f,"snp.id", val=snp,compress="ZIP",closezip=TRUE)
add.gdsn(dat.f,"snp.chromosome", val=chr,compress="ZIP",closezip=TRUE)
add.gdsn(dat.f,"snp.position", val=pos,compress="ZIP",closezip=TRUE)
add.gdsn(dat.f,"snp.allele", val=all,compress="ZIP",closezip=TRUE)
add.gdsn(dat.f,"genotype", val=gen,compress="ZIP",closezip=TRUE,storage="bit2",replace=TRUE)
add.gdsn(dat.f,"sample.annot", val=ant,compress="ZIP",closezip=TRUE)
put.attr.gdsn(index.gdsn(dat.f,"genotype"),"sample.order")
closefn.gds(dat.f)	
# add coordinates
geo<-read.table("Geo_SPGR.csv",header=TRUE,sep=";",dec=",")


