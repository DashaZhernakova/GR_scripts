library(pophelper)
alist <- readQ('OneDrive/Python-practice/admixture-GR/fltData.4.Q')

twolabset <- read.delim('OneDrive/Python-practice/admixture-GR/popRegion.txt', header=T,stringsAsFactors=F)

onelabset <- twolabset[,1,drop=FALSE]

plotQ(alist, grplab=twolabset, 
             selgrp="Region", 
             outputfilename='admixPlot',
             sortind="all",
             grplabangle=-90, 
             grplabjust=0.5,
             grplabheight=3,
             height=2,
             width=25,)