library(ggplot2)
d <- read.table("/Users/dashazhernakova/Documents/Doby/Bake-off/report_november_2016/HLA_typing/HLA_platform_cmp2.txt", header = T, as.is = T, check.names = F)
d$gene = as.factor(d$gene)
d$platform = as.factor(d$platform)

# workaround to set the limits for vertical grey lines (geom_linerange)
for (g in unique(d$gene)){ 
  d[d$gene == g,"lower"] <- min(d[d$gene == g, "fraction_errors"])
  d[d$gene == g,"upper"] <- max(d[d$gene == g, "fraction_errors"])
}

pdf("/Users/dashazhernakova/Documents/Doby/Bake-off/report_november_2016/HLA_typing/HLA_platform_cmp2.pdf", useDingbats = F)
ggplot() +
  geom_linerange(data=d, mapping=aes( x = gene, ymin = lower, ymax = upper), color = "grey", size = 1.7) +
  geom_point(data=d, mapping=aes(x = gene, y = fraction_errors), size = 3) +
  xlab("gene") +
  scale_y_continuous(name = "Percent errors", limits = c(0,1), labels = scales::percent_format(accuracy = 1)) + 
  theme_bw()
dev.off()
