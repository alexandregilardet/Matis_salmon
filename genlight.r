library(adegenet)

# import
g <- fasta2genlight(
'C:/Users/alexa/Matis/data/SSaTrack_Island_Stofnfiskur_fullMT.fa'
, chunkSize = 10, parallel = FALSE)
g <- read.PLINK(
        file = "C:/Users/alexa/Matis/results/20_04_20/ssa09_recodeA_plink/SSaTrack_Island_chr09_recodeA.raw",
        map.file = "C:/Users/alexa/Matis/results/20_04_20/ssa09_recodeA_plink/SSaTrack_Island_chr09.map",
        quiet = FALSE,
        chunkSize = 1000,
        parallel = FALSE,
        )

locNames(g)# position.alleles

# ssa25 51481326
akap11 <- 28712315:28742294
vgll3 <- 28654947:28659019
akap11vgll3 <- 28654947:28742294

# ssa09 141712163
six6 <- 24902777:24905552

# plot position in aln
temp <- density(position(g), bw=100000)
plot(temp, type="n", xlab="Position in the alignment",
     main="Location of the SNPs", xlim=c(0,141712163))
polygon(c(temp$x,rev(temp$x)), c(temp$y, rep(0,length(temp$x))),
        col=transp("blue",.3))
points(position(g), rep(0, nLoc(g)), pch="|", col="blue")
# work bandwidth, smaller is sharper
points(six6,rep(0, length(six6)) , pch="|", col="red")

# plot position in genome
snpposi.plot(position(g), genome.size=141712163, codon=FALSE)
#test polymorphic hotspot
snpposi.test(position(g), genome.size=141712163)

# import saving allele nb at all pos in other slot
#g2 <- fasta2genlight('C:/Users/alexa/Matis/data/SSaTrack_Island_Stofnfiskur_fullMT.fa',
       #chunk=10,saveNbAlleles=TRUE, quiet=TRUE, parallel=FALSE)
# percentage polymorphic sites
100*mean(unlist(other(g2))>1)
# allele nb distribution
#temp <- table(unlist(other(g2)))
#barplot(temp, main="Distribution of the number nnof alleles per loci",
        #xlab="Number of alleles", ylab="Number of sites", col=heat.colors(4))

glPlot(g, posi="topleft")
# white missing data

myFreq <- glMean(g)
myFreq <- c(myFreq, 1-myFreq)
hist(myFreq, proba=TRUE, col="darkseagreen3", xlab="Allele frequencies",
     main="Distribution of allele frequencies", nclass=20)
temp <- density(myFreq, bw=.05)
lines(temp$x, temp$y*2,lwd=3)
# fixed loci 0 and 1, intermediate frequencies

temp <- density(glNA(g), bw=10)
plot(temp, type="n", xlab="Position in the alignment", main="Location of the missing values", xlim=c(0,16665))
polygon(c(temp$x,rev(temp$x)), c(temp$y, rep(0,length(temp$x))), col=transp("blue",.3))
points(glNA(g), rep(0, nLoc(g)), pch="|", col="blue")
# beginning probably reflecting heterogeneity in DNA amplication during sequencing

# PCA
pca1 <- glPca(g)
scatter(pca1, posi="bottomright")
title("PCA of the data\n axes 1-2")