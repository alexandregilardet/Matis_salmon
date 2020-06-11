library(plyr)

# Prepare input file for the mitochondrial reference genome to use in final_SNP_on_mt.R
# Only main genes are conserved for the plot 

read.csv("C:/Users/alexa/Matis/results/26_02_20/NC_001960_gene_lengths.csv") -> mtref

gsub("gene", "", mtref$Name) -> mtref$Name
rename(mtref, c("Name" = "Gene")) -> mtref
as.numeric(gsub(",", "", mtref$Length)) -> mtref$Length

mtref[!grepl('CDS', mtref$Gene),] -> mtref

mtref[!grepl('Pro', mtref$Gene),] -> mtref
mtref[!grepl('Thr', mtref$Gene),] -> mtref
mtref[!grepl('Leu', mtref$Gene),] -> mtref
mtref[!grepl('Ser', mtref$Gene),] -> mtref
mtref[!grepl('His', mtref$Gene),] -> mtref
mtref[!grepl('Arg', mtref$Gene),] -> mtref
mtref[!grepl('Lys', mtref$Gene),] -> mtref
mtref[!grepl('Asp', mtref$Gene),] -> mtref
mtref[!grepl('Tyr', mtref$Gene),] -> mtref
mtref[!grepl('Cys', mtref$Gene),] -> mtref
mtref[!grepl('Asn', mtref$Gene),] -> mtref
mtref[!grepl('Ala', mtref$Gene),] -> mtref
mtref[!grepl('Trp', mtref$Gene),] -> mtref
mtref[!grepl('Met', mtref$Gene),] -> mtref
mtref[!grepl('Gln', mtref$Gene),] -> mtref
mtref[!grepl('Val', mtref$Gene),] -> mtref
mtref[!grepl('Phe', mtref$Gene),] -> mtref

mtref[!grepl('origin', mtref$Gene),] -> mtref

mtref[c("Length", "Gene")] -> mtref # invert columns

1:nrow(mtref) -> row.names(mtref)# reinitialize row numbers after removal of some

write.csv(mtref,"C:/Users/alexa/Matis/results/26_02_20/mtref.csv")