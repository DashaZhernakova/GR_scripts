args <- commandArgs(trailingOnly = TRUE)

library(gdsfmt)
library(SNPRelate)

path.fn <- args[1]
gds_data <- openfn.gds(paste(path.fn,".gds",sep=""))
snpgdsGDS2BED(gds_data, path.fn)

