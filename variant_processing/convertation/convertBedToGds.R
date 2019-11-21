args <- commandArgs(trailingOnly = TRUE)

library(gdsfmt)
library(SNPRelate)

path.fn <- args[1]

snpgdsBED2GDS(paste(path.fn,".bed",sep=""),paste(path.fn,".fam",sep=""),paste(path.fn,".bim",sep=""),paste(path.fn,".gds",sep=""))

