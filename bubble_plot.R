library(ggplot2)
library(plyr)
library(dplyr)

fishing <- fishing_counts_2016
gts <- rename(gts, c("River_short" = "River"))
data <- merge(fishing, gts, by = c("River_short"))
data$River_short <- gsub("LaD", "La_D", data$River_short)

tbl_df = data.frame()
#for global plot
tbl_list = list()
write.table(tbl_list, file = "near_vgll3_iHS_bubble.csv", 
            col.names = NA, row.names = TRUE, sep = ",")
#to have empty column for row name
plot_list = list() 
#empty

for (i in 1:nrow(data)) {
  row <- data[i,] #row by row
  row.names(row) <- as.character(row$River_short) 
  #remove levels attribute
  
  river <- as.character(row$River_short) 
  SW <- row$freq_2SW
  nd_1 <- colnames(row)[6]
  nd_2 <- colnames(row)[7]
  nd_3 <- colnames(row)[8]
  gt_1 <- pull(row[6])
  gt_2 <- pull(row[7])
  gt_3 <- pull(row[8])
  nd <- c(nd_1, nd_2, nd_3)
  gt <- c(gt_1, gt_2, gt_3)
  freq_2SW <- c(SW, SW, SW)
  
  #create output table
  tbl <- select(.data = row, freq_2SW, 
                       nd_1, nd_2, nd_3)
  tbl_list[[i]] = tbl
  write.table(tbl, file = "near_vgll3_iHS_bubble.csv", 
              append = TRUE,col.names = NA, row.names = TRUE,
              sep = ",")
  
  #create table for the plot
  input <- data.frame(nd, freq_2SW, gt)
  
  #global table for global plot
  tbl_df <- rbind(input, tbl_df)
  
  #bubble plot
  plot <- ggplot(input, 
    aes(x = nd, y = freq_2SW, size = gt))+    
    geom_point()+
    ylim(0, 0.70)+
    scale_size_continuous(limits=c(0,1),
          breaks=c(0.05,0.25,0.5,0.75,1))+
    theme_minimal()+
    ggtitle(river)
  plot_list[[i]] = plot
  
  file_name = paste("bubble_plot_", river, ".png", sep="")
  png(file_name)
  print(plot_list[[i]])
  dev.off()
}


glob_plot <- ggplot(tbl_df, aes(x = nd, 
                               y = freq_2SW, 
                               size = gt))+    
  geom_point(alpha = 0.3)+
  ylim(0, 0.70)+
  scale_size_continuous(limits=c(0,1),
                        breaks=c(0.05,0.25,0.5,0.75,1))+
  theme_minimal()

png("bubble_plot_global.png")
print(glob_plot)
dev.off()

#put SNP of interest into gt
#choose the name for csv file twice
#replace 0 by 0.00 in csv table