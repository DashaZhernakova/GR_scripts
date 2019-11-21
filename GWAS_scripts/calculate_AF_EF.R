library(SNPRelate)
#setwd("/Users/dashazhernakova/Documents/Doby/HIV/")

args <- commandArgs(trailingOnly = TRUE)

# calculates attiributable fraction and explained fraction according to Nelson and O'Brien 2006
calculate_AR_EF <- function(tab){
  #relative risk
  RR <- (tab[2,1]/sum(tab[2,])) / (tab[1,1]/sum(tab[1,]))
  # frequency of exposure
  f <- sum(tab[2,])/(sum(tab[1,]) + sum(tab[2,]))
  # Attributable fraction
  AR <- f*(RR-1)/(1+f*(RR-1))
  
  # Entropy
  #E_exposure <- -1*( (sum(tab[1,])/sum(tab)) * log(sum(tab[1,])/sum(tab)) + (sum(tab[2,])/sum(tab)) * log(sum(tab[2,])/sum(tab)))
  #E_disease <- -1*( (sum(tab[,1])/sum(tab)) * log(sum(tab[,1])/sum(tab)) + (sum(tab[,2])/sum(tab)) * log(sum(tab[,2])/sum(tab)))
  #E_joint <- -1*( (sum(tab[1,1])/sum(tab)) * log(sum(tab[1,1])/sum(tab)) + (sum(tab[2,1])/sum(tab)) * log(sum(tab[2,1])/sum(tab)) + (sum(tab[1,2])/sum(tab)) * log(sum(tab[1,2])/sum(tab)) + (sum(tab[2,2])/sum(tab)) * log(sum(tab[2,2])/sum(tab)))
  
  #E_joint <- Entropy(tab,base=exp(1))
  
  # Mutual information
  marg_i <- c(sum(tab[1,])/sum(tab), sum(tab[2,])/sum(tab))
  marg_j <- c(sum(tab[,1])/sum(tab), sum(tab[,2])/sum(tab))
  
  freq_tab <- tab/sum(tab)
  
  I <- freq_tab[1,1]*log(freq_tab[1,1]/(marg_i[1]*marg_j[1])) + freq_tab[1,2]*log(freq_tab[1,2]/(marg_i[1]*marg_j[2])) + freq_tab[2,1]*log(freq_tab[2,1]/(marg_i[2]*marg_j[1])) + freq_tab[2,2]*log(freq_tab[2,2]/(marg_i[2]*marg_j[2])) 
  # I <- MutInf(freq_tab)
  
  Imax <- -1*(marg_j[1]*log(marg_j[1]) + marg_j[2]*log(marg_j[2]))
  
  EF <- I/Imax
  #EF
  return(c(RR, AR, EF))
}

#fname="array.original.only.nomale.maf0.05.hwe0.0001.geno0.05.gds"
fname = args[1]
dat.f<-openfn.gds(fname,readonly=FALSE)
gen<-read.gdsn(index.gdsn(dat.f,"genotype"))
id<-read.gdsn(index.gdsn(dat.f,"sample.id"))
snp<-read.gdsn(index.gdsn(dat.f,"snp.id"))
chr<-read.gdsn(index.gdsn(dat.f,"snp.chromosome"))
pos<-read.gdsn(index.gdsn(dat.f,"snp.position"))
all<-read.gdsn(index.gdsn(dat.f,"snp.allele"))
ant<-read.gdsn(index.gdsn(dat.f,"sample.annot"))
closefn.gds(dat.f)

hit_table_fname = args[2]
#hit_table_fname = "/Users/dashazhernakova/Documents/Doby/HIV/MA_hits.txt"
hit_table <- read.table(hit_table_fname, header = T, as.is = T, check.names = F)

hit_table["RR"] <- NA
hit_table["AR"] <- NA
hit_table["EF"] <- NA


for (hit_num in 1:nrow(hit_table)){
  chr_pos <- paste0(hit_table[hit_num,"CHR"], ":", hit_table[hit_num,"BP"])
  #chr_pos = "17:62006259"
  #print(chr_pos)
  
  snp_idx <- which(snp==chr_pos)
  if (length(snp_idx) < 1){
    print("No such SNP! SNP are named in the following way in this dataset:")
    paste(head(snp))
  }
  
  geno <- as.numeric(gen[, snp_idx])
  flt <- geno != 3 # remove missing genotypes
  pht<-(ant$phenotype-1)[flt]
  
  geno_table0 <- as.matrix(table(geno[flt],pht))
  
  # fill empty genotype counts with 0
  if (! "0" %in% rownames(geno_table0)){
    geno_table0 <- rbind(c(0,0), geno_table0)
    row.names(geno_table0)[1] <- "0"
  }
  if (! "1" %in% rownames(geno_table0)){
    geno_table0 <- rbind(geno_table0, c(0,0))
    row.names(geno_table0)[3] <- "1"
  }
  if (! "2" %in% rownames(geno_table0)){
    geno_table0 <- rbind(geno_table0, c(0,0))
    row.names(geno_table0)[3] <- "2"
  }
  
  geno_table <- matrix(nrow = 3, ncol = 2)
  colnames(geno_table) <- c("case", "control")
  row.names(geno_table) <- row.names(geno_table0)
  #print(geno_table0)
  geno_table[,1] <- geno_table0[,2]
  geno_table[,2] <- geno_table0[,1]
  #geno_table
  
  dom_table <- matrix(nrow = 2, ncol = 2)
  dom_table[1,] <- geno_table["0",] + geno_table["1",]
  dom_table[2,] <- geno_table["2",]
  #dom_table
  
  rec_table <- matrix(nrow = 2, ncol = 2)
  rec_table[1,] <- geno_table["0",]
  rec_table[2,] <- geno_table["1",] + geno_table["2",]
  #rec_table
  
  allelic_table <- matrix(nrow = 2, ncol = 2)
  allelic_table[1,] <- 2*geno_table["0",] + geno_table["1",]
  allelic_table[2,] <- geno_table["1",] + 2*geno_table["2",]
  #allelic_table
  
  #res_table <- matrix(nrow = 3, ncol = 3)
  #colnames(res_table) <- c("RR", "AF", "EF")
  #row.names(res_table) <- c("DOM", "REC", "ALL")
  #res_table[1,] <- calculate_AR_EF(dom_table)
  #res_table[2,] <- calculate_AR_EF(rec_table)
  #res_table[3,] <- calculate_AR_EF(allelic_table)
  #print(hit_num)
  if (hit_table[hit_num,"TEST"] == "DOM"){
    hit_table[hit_num,c("RR", "AR", "EF")] <- calculate_AR_EF(dom_table)
  } else if (hit_table[hit_num,"TEST"] == "REC"){
    hit_table[hit_num,c("RR", "AR", "EF")] <- calculate_AR_EF(rec_table)
  } else if (hit_table[hit_num,"TEST"] == "ALLELIC"){
    hit_table[hit_num,c("RR", "AR", "EF")] <- calculate_AR_EF(allelic_table)
  } else {
    hit_table[hit_num,c("RR", "AR", "EF")] <- "NA"
  }
  
}

write.table(hit_table, file = paste0(hit_table_fname, ".AR_EF.txt"), sep = "\t", quote = F,  row.names = F)




