library(pophelper)

# Variables

title_label <- 'Sorted Ancient and Modern'
admx_path <- '~/Downloads/tmp/sorted_GR+papers.hg19+ancient.shared_snps.Q'
fam_path <- '~/Downloads/tmp/sorted_GR+papers.hg19+ancient.shared_snps.txt'
out_path <- '~/Downloads/tmp/sorted_GR+papers.hg19+ancient.shared_snps'

# Plotting

admx_results <- readQ(admx_path)

sample_info <- read.delim(fam_path,
                          sep = ' ',
                          header = F,
                          stringsAsFactors = F)

family_labels <- sample_info[, 1, drop = F]
sample_labels <- sample_info$V2

rownames(admx_results[[1]]) <- sample_labels

plotQMultiline(admx_results,
               lpp = 5,
               showtitle = T,
               titlelab = title_label,
               showsubtitle = T,
               showindlab = T,
               useindlab = T,
               grplab = family_labels,
               outputfilename = out_path
)
