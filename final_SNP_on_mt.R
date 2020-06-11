
library(plyr)
library(BioCircos)

# Plot fixed mitochondrial SNPs from both populations

# Prepare genome

read.csv("C:/Users/alexa/Matis/results/26_02_20/mtref.csv") -> mtref

NULL -> mtref$X # remove 1st column with numbers of row after save and re-open

split(mtref$Length, mtref$Gene) -> mtref # into list

# Prepare total mitochondrial SNPs to get the gene and coordinate inside
# according to position of later selected SNPs

read.csv("C:/Users/alexa/Matis/results/26_02_20/total_snp_coord.csv", sep=";") -> snp

rename(snp, c("�..snp_gene" = "snp_gene")) -> snp
gsub("gene", "", snp$snp_gene) -> snp$snp_gene

# Prepare combined SNPs
read.csv("C:/Users/alexa/Matis/results/25_03_20/SSaTrack_fullMT_MAF_combined.csv", header = TRUE, sep = ";") -> maf

rename(maf, c("Position" = "pos")) -> maf
gsub(",", ".", maf$IS.MAF) -> maf$IS.MAF
gsub(",", ".", maf$Stofn.MAF) -> maf$Stofn.MAF

merge(snp, maf) -> snp_maf

as.numeric(snp_maf$snp_coord) -> snp_maf_coord
as.character(snp_maf$snp_gene) -> snp_maf_gene
as.numeric(snp_maf$IS.MAF) -> values_i_maf
as.numeric(snp_maf$Stofn.MAF) -> values_s_maf

# Scale 

min_val = 0
max_val = 1

# Add fixed SNPs

read.csv("C:/Users/alexa/Matis/results/25_03_20/fixed_snp.csv", header = TRUE, sep = ";") -> fixed

merge(snp, fixed) -> snp_fixed

as.numeric(snp_fixed$snp_coord) -> fixed_coord
as.character(snp_fixed$snp_gene) -> fixed_gene
as.numeric(snp_fixed$IS.MAF) -> fixed_i_maf
as.numeric(snp_fixed$Stofn.MAF) -> fixed_s_maf

# Plot

tracks = BioCircosSNPTrack('SNP_Island', snp_maf_gene, snp_maf_coord, values_i_maf, size = 1.5, color = 'red', opacities = 1, range = c(min_val, max_val))
tracks = tracks + BioCircosSNPTrack('SNP_Stofnfiskur', snp_maf_gene, snp_maf_coord, values_s_maf, size = 1.5, color = 'darkblue', opacities = 0.6, range = c(min_val, max_val))
tracks = tracks + BioCircosSNPTrack('Fixed_Island', fixed_gene, fixed_coord, fixed_i_maf, size = 5, color = 'red', opacities = 0.5, range = c(min_val, max_val))
tracks = tracks + BioCircosSNPTrack('Fixed_Stofnfiskur', fixed_gene, fixed_coord, fixed_s_maf, size = 5, color = 'darkblue', opacities = 0.5, range = c(min_val, max_val))

tracks = tracks + BioCircosBackgroundTrack("myBackgroundTrack", borderColors = "#AAAAAA", borderSize = 0.6)

BioCircos(tracks, genome = mtref, displayGenomeBorder = FALSE, genomeFillColor = c("red", "darkblue"), zoom = TRUE, chrPad = 0.05, genomeTicksDisplay = FALSE, genomeLabelTextSize = 9, genomeLabelOrientation = 90, genomeLabelDy = 35, height = 500)


