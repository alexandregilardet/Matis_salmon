data <- cor_per_river

Bla <- data[1:5]
Bre <- data[6:10]
Ell <- data[11:15]
Gri <- data[16:20]
Haf <- data[21:25]
Hof <- data[26:30]
Hvi <- data[31:35]
Kro <- data[36:40]
LaA <- data[41:45]
LaD <- data[46:50]
Lad <- data[51:55]
Lan <- data[56:60]
Lau <- data[61:65]
Nor <- data[66:70]
Olf <- data[71:75]
Sun <- data[76:80]
Thj <- data[81:85]
Vid <- data[86:90]

cc = cor(data, method = 'pearson')

# heatmap corrplot
#library(corrplot)
#corrplot(cc)
# get clusters of similar patterns
#corrplot(cc, tl.col = "black", order = "hclust", 
         #hclust.method = "average", addrect = 4, 
         #tl.cex = 0.7)
# save
#png(filename = "Bla_corrplot.png", width = 6, height = 6, 
    #units = "in", res = 400)

# heatmap ggplot
#library(ggplot2)
#library(reshape2)
cc_df = as.data.frame(cc)
cc_df$snp = row.names(cc_df)
ccm = melt(cc_df, id = "snp") # long format
# keep same initial order
ccm$snp <- factor(ccm$snp,levels=unique(ccm$snp))
xx = ggplot(ccm, aes(x = variable, y = snp)) + 
  geom_tile(aes(fill = value), colour = "grey45") + 
  coord_equal() + 
  scale_fill_gradient(low = "navy", high = "red") + 
  theme(axis.text.y = element_text(size =2, face = "bold", colour = "grey25"), 
        legend.title = element_text(size = 10, face = "bold"),legend.position = "bottom", 
        axis.text.x = element_text(size = 2, angle = 90, face = "bold",colour = "grey25", vjust = 0.5, hjust = 0), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = NA), 
        axis.ticks = element_blank()) + 
  labs(x= "", y = "", fill = "Pearson's Correlation") + 
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev(levels(ccm$snp))) 
xx
# save
ggsave("Total_heatmap.png")
dev.off()
