library("reshape2")
setwd("/Users/dashazhernakova/Documents/Doby/GenomeRussia/freeze2019/population_genetics/ibd/")
d <- read.table("all_chr.GR+papers.eurasian.subset.ibd.format.segments.step0.5", header = T, as.is = T, check.names = F)
d2 <- melt(d, id.vars = c("Population1", "Population2"), measure.vars =colnames(d)[5:ncol(d)])
#pop_subset <- c("Forest-Nenets","Komis","Maris","Evenks","Voronezh","Bashkirs","Vepsas","Tundra-Nenets","Khantys","Mansis","Finnish","Ingrians","Estonians","Saami","Arkhangelsk","Novgorod","Pskov","Tatars","Udmurts","Mordvins","Chuvashes")
pop_subset <- c("Karelians","Vepsas", "Komis","Finnish","Latvians","Lithuanians","Estonians","Swedes","Voronezh", "Udmurts")
pop_subset <- c("Karelians","Komis","Finnish","Vepsas","Ingrians", "Saami", "Latvians","Lithuanians","Estonians","Voronezh", "Udmurts", "Mansis", "Pskov", "Novgorod")
d2_subs <- d2[d2$Population2 %in% pop_subset, ]

gr_pops = c("Pskov","Novgorod", "Arkhangelsk","Rus_Udmurtia","Vladimir","Yaroslavl","Rostov","Voronezh")
plots <- list()

cnt = 1
for (pop in pop_subset){
print (pop)
  pop_d <- d2_subs[d2_subs$Population1 == pop,]
  p <- ggplot(data = pop_d, aes(x=variable, y = value, color = Population2, group = Population2)) +
  geom_point(aes(shape = Population2)) +
  geom_line(size=0.5) +
  ggtitle(pop) + 
  labs(x = "IBD seg length", y = "# segments") +
  scale_shape_manual(values=seq(1,16)) +
    ylim(0,max(6,max(pop_d$value))) +
    theme(plot.title = element_text(size=20))
  plots <- c(plots, list(p))
 cnt = cnt + 1
}

pdf("ibd_fu.pdf", height = 8, width = 15)
multiplot(plots)
dev.off()

pop_subset <- c("Komis","Voronezh","Vepsas","Mansis","Finnish","Ingrians","Estonians","Saami","Arkhangelsk","Udmurts", "Karelians")
pop_subset <- c("Komis","Voronezh","Vepsas","Mansis","Finnish","Ingrians","Estonians","Saami","Arkhangelsk","Udmurts", "Karelians", "Khantys", "Pskov", "Novgorod", "Yaroslavl", "Swedes")
#pop_subset <- c("Komis","Voronezh","Bashkirs","Tatars","Udmurts","Mordvins","Chuvashes")

pop_subset <- c("Komis","Vepsas","Mansis","Finnish","Ingrians","Estonians","Saami","Udmurts", "Karelians", "Khantys", "Latvians")
pop="Estonians"

pdf(paste0(pop,".pdf"), height = 8, width = 15)
d2_subs <- d2[d2$Population2 %in% pop_subset, ]
pop_d <- d2_subs[d2_subs$Population1 == pop,]

ggplot(data = pop_d, aes(x=variable, y = value, color = Population2, group = Population2)) +
  geom_point(aes(shape = Population2)) +
  geom_line(size=0.5) +
  ggtitle(pop) + 
  labs(x = "IBD seg length", y = "# segments") +
  scale_shape_manual(values=seq(1,22)) +
  ylim(0,6) +
  theme(plot.title = element_text(size=20))
dev.off()

pdf("ibd_rus-udm.pdf", height = 8, width = 15)
p2
dev.off()

col4 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F",
                           "cyan", "#007FFF", "blue", "#00007F"))

pdf("corrplots_ibd3.pdf", height = 10, width = 15, useDingbats=F)
par(mfrow=c(4,2))
for (pop in gr_pops){
  print (pop)
pop_d <- d2_subs[(d2_subs$Population1 == pop),]
pop_d_matrix <- acast(pop_d, Population2~variable, value.var="value")

pop_d_matrix[which(pop_d_matrix > 6)] = 6
pop_d_matrix[which(pop_d_matrix < -6)] = -6
corrplot(as.matrix(pop_d_matrix), tl.cex=1, tl.col=gray(0.3), type='full', 
         cl.lim=c(0, 6), is.corr=F, rect.lwd=0.2,  tl.offset = 0.5, 
         cl.cex = 0.8, cl.ratio = 0.3, main = pop, col = col4(100), mar=c(0,0,1,0))

}
dev.off()






multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
