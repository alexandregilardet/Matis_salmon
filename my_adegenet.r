library(adegenet)
library(ape)
library(hierfstat)
library(pegas)
library(poppr)

#Imports
names(d) <- gsub("\\.", "_", names(d)) #remove dots
x<-df2genind(d[4:38], ploidy=2, ind.names=d$short_name, sep = "", NA.char="NN")
#x<-df2genind(d[3:126], ploidy=2, ind.names=d$ind, pop=d$pop, sep=":", NA.char="NA")
#x <- read.genetix(file="C:/Users/alexa/Matis/data/David_17_03_20/AlexPOP2GTX-GENETIX.gtx")
x<-df2genind(d[3:11], ploidy=2, ind.names=d$ind, pop=d$pop, sep="", NA.char="N")

#Presence/absence @type not codom but PA
d <- data.frame(lapply(d, function(x) 
{gsub("NN", "0", x)})) #NN into 0
d[4:21] <- data.frame(lapply(d[4:21], function(x) 
{as.integer(x!="0")}))
x<-df2genind(d[4:5], ploidy=1, ind.names=d$short_name, type = "PA")

pca1 <- dudi.pca(x, scale = FALSE, scannf = FALSE)
plot(pca1$li, cex = 1)
abline(h=0,v=0,col="grey",lty=2)
s.arrow(pca1$c1*.2, add.plot=TRUE)
     
#Basics
informloci(x) -> x #uninformative loci, drop them
summary(x) -> toto 
names(toto)

#x <- x[,loc=c('AX_87532311', 'AX_96317592', 'AX_87665320', 'AX_87220302', 'AX_87268437', 'AX_87570456')] #select loci

plot(toto$n.by.pop, toto$pop.n.all, xlab="Colonies sample size",
     ylab="Number of alleles",main="Alleles numbers and sample sizes",
     type="n")
text(toto$n.by.pop,toto$pop.n.all,lab=names(toto$n.by.pop))

barplot(toto$loc.n.all, ylab="Number of alleles",
        main="Number of alleles per locus")

barplot(toto$Hexp-toto$Hobs, main="Heterozygosity: expected-observed",
        ylab="Hexp - Hobs")

barplot(toto$n.by.pop, main="Sample sizes per population",
        ylab="Number of genotypes",las=3)

#Is mean observed H significantly lower than mean expected H?
bartlett.test(list(toto$Hexp,toto$Hobs))

#reject homogeneity if p < .05
t.test(toto$Hexp,toto$Hobs,pair=T,var.equal=TRUE,alter="greater")
#significative difference if p < .05

#HWT
x.hwt <- hw.test(x, B=0) #0 for parametric otherwise second test Monte-Carlos 
x.hwt #each locus
#reject equilibrium if Pr < .05

#F statistics
fstat(x) #overall
#Fst (pop/total), Fit (Ind/total), Fis(ind/pop)
Fst(as.loci(x))

#Are these values significant?
Gtest <- gstat.randtest(x, nsim = 99)
Gtest
plot(Gtest)

matFst <- pairwise.fst(x)
matFst
is.euclid(matFst)
#PCA requires Euclideanity

#PCA
sum(is.na(x$tab)) #nb NA
X <- tab(x, freq = TRUE, NA.method = "mean")
#tab is allele counts/frequencies
#replace NAs with mean allele frequency
pca1 <- dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 3)
#disable scaling bc all alleles vary on same scale
#remove two last args for interactive selection nb axes
barplot(pca1$eig[1:50], main = "PCA eigenvalues", col = heat.colors(50))
names (pca1)
#principal components of analysis ; synthetic variables summarizing genetic diversity
s.label(pca1$li)
title("PCA of X dataset\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

s.class(pca1$li, pop(x))
title("PCA of X dataset\naxes 1-3")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

#distance and color represent genetic diff
colorplot(pca1$li, pca1$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 2")
title("PCA of x dataset\naxes 1-2")
abline(v=0,h=0,col="grey", lty=2)

colorplot(pca1$li[c(1,3)], pca1$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 3")
title("PCA of x dataset\naxes 1-3")
abline(v=0,h=0,col="grey", lty=2)

#$c1 allele loadings

#DAPC
#find clusters
grp <- find.clusters(x)
#keep maximum PCs to keep all information 
#elbow matches nb clusters 
names(grp)
table(pop(x), grp$grp) #access actual groups with pop
table.value(table(pop(x), grp$grp), col.lab=paste("inf", 1:2),
            row.lab=paste("ori", 1:2))

#browser version but first change data to genetix format
adegenetServer(what=c("DAPC"))

#manual
#uses pre-defined grps!
dapc1 <- dapc(x)
#retain less PCs is better 
#ok to retain all discriminants if few
scatter (dapc1)
scatter (dapc1, posi.da="bottomright", scree.pca=TRUE,
         posi.pca="bottomleft")

#contribution of alleles
contrib <- loadingplot(dapc1$var.contr, threshold = 0.015)

#grp membership #prior/posterior
assignplot(dapc1, subset =) 
#blue cross is prior pre-defined attribution
