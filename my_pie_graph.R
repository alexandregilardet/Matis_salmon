#Making pie charts from table, works with NAs

library(ggplot2)
library(reshape)
library(reshape2)
library(ggsci)
library(RColorBrewer)

#colour blind-friendly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#select snp concerned in total data
d_tot <- data[c(1,81:84)]
d_tot
d <- d_tot[1:(length(d_tot)-1)] #drop last col tot

id <- colnames(d)[2]
marker <- strsplit(id, "_")[[1]][1] 
marker
#subset list then vector in list

#new col fusing river and number of ind
d$River <- paste(d_tot$River,d_tot[,ncol(d_tot)], sep="")
d

#Table is river, gt_freq1, gt_freq2, gt_freq3

d2 <- melt(d) #change the table to long form
d2 <- rename(d2, c("variable" = "Genotype"))

colourCount = length(unique(d2$Genotype)) #to count the colours to use
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

#River is id variables
q1 <- ggplot(d2, aes(x="", y=value, group=Genotype, color=Genotype, fill=Genotype)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) + facet_wrap(~ River) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank())
  #ggtitle(marker)

#q2 <- q1 + scale_fill_manual(values = getPalette(colourCount))
q2 <- q1 + scale_fill_manual(values = cbPalette)


fname = paste("pie_", marker, ".png", sep="")
png(fname)
print(q2)
dev.off()