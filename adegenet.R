##For all populations

library(adegenet)
library(dplyr)
library(tidyr)
library("hierfstat")
library(spdep)
library(readxl)
library(adespatial)
library(ggplot2)
library(RColorBrewer)
library(ade4)

# Load csv file
microsat_data <- read.csv("DatasetR_2.csv", header = TRUE)
microsat_data

# number of individuals
nind=nrow(microsat_data)
# number of loci, by removing columns of info
nloci=ncol(microsat_data)-7

# Convert to genind object using the combined data and specifying the separator
genind_obj <- df2genind(microsat_data[,c(8:14)], sep = "/", NA.char = "NA")
# Inspect the genind object
summary(genind_obj)

# add the (sub)population information
pop(genind_obj) <- microsat_data[,3]
genind_obj

# add geographic coordinates and create a first map
genind_obj@other$xy <-microsat_data[,5:6]
print(genind_obj@other$xy)

png("Geographic coordinates.png", width = 8, height = 6, units="in", res=1080)
plot(genind_obj@other$xy, cex=1.5, xlab='x', ylab='y')
points(genind_obj@other$xy, col="red",pch=20)

dev.off()


# Calculate some basic genetic parameters for the whole dataset
sum_obj <-summary(genind_obj)
sum_obj

# Expected and observed hererozygosity and allelic richness
Ho_all <- mean(sum_obj$Hobs)
He_all <- mean(sum_obj$Hexp)
F_all <- (He_all-Ho_all)/He_all
Na_all <- mean(sum_obj$loc.n.all)
Na_all
Ho_all
He_all
F_all

summary_table <- data.frame(
  Na_all,
  Ho_all,
  He_all,
  F_all
)

print(summary_table)

##write.table(summary_table, "Expected and observed hererozygosity and allelic richness.txt", sep = "\t", row.names = FALSE, quote = FALSE)


# Convert genind object to genpop object
genpop_obj <- genind2genpop(genind_obj)
genpop_obj

popNames(genpop_obj)

summary(genpop_obj)


## Plot - Sample size per Population

png("Sample Sizes per Pop.png", width = 7, height = 5, units="in", res=1080)
colors <- brewer.pal(n = 7, name = "Paired")

# Create the barplot with customized aesthetics
barplot(sum_obj$n.by.pop,
        main = "Sample Sizes per Population",
        ylab = "Number of Genotypes",
        xlab = "Population",
        col = colors,               # Set bar colors
        border = "black",           # Set bar borders to black
        las = 1,                    # Rotate x-axis labels horizontally
        cex.names = 0.8,            # Reduce font size of x-axis labels
        cex.lab = 1.2,              # Increase font size of axis labels
        cex.main = 1.4,             # Increase font size of title
        # Expand y-axis limits slightly
        ylim = c(0, max(sum_obj$n.by.pop) * 1.1)  )
dev.off()


# F-statistics: Fst and Fis
wc(genind_obj)

FST <- wc(genind_obj)$FST
FIS <- wc(genind_obj)$FIS 

wc_table <- data.frame(
  FST,
  FIS
)

print(wc_table)

##write.table(wc_table, "Fst and Fis.txt", sep="\t", row.names = FALSE, quote = FALSE)

# Genetic distances (Nei) between populations
nei <- genet.dist(genind_obj, method = "Nei87")
nei1 <- genet.dist(genind_obj, method = "Dm")
nei1

# Convert the 'nei' object to a data frame (if it's a matrix or dist object)
nei_df <- as.data.frame(as.matrix(nei1))
nei_df


# Perform UPGMA hierarchical clustering
hc_upgma <- hclust(as.dist(nei1), method = "average")

# Plot the UPGMA dendrogram
png("UPGMA Dendrogram_all.png", width = 6 , height= 4.5, units= "in", res = 1080)
plot(hc_upgma, main = "UPGMA Dendrograms of Genetic Distances (Nei's Dm)", xlab = "Samples", ylab = "Genetic Distance")

dev.off()

# Optional: 'ape' to create a more advanced tree plot
phylo_tree <- as.phylo(hc_upgma)
plot(phylo_tree, type = "fan", main = "UPGMA Dendrogram (Circular Layout)")



# Save the Nei's distance table to a .txt file
###write.table(nei_df, "nei_genetic_distance.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)


# New genind object out of obj, where populations are separated
newobj<-seppop(genind_obj)
names(newobj)

obj1<-newobj$Dafni
obj2<-newobj$Panteleimon
obj3<-newobj$Center
obj4<-newobj$Vela
obj5<-newobj$Psaria
obj6<-newobj$Perama
obj7<-newobj$Maries
sum1<-summary(obj1)
sum2<-summary(obj2)
sum3<-summary(obj3)
sum4<-summary(obj4)
sum5<-summary(obj5)
sum6<-summary(obj6)
sum7<-summary(obj7)

table_Heterozygosity <- data.frame(
  Population = c("Dafni", "Panteleimon", "Center", "Vela", "Psaria", "Perama", "Maries"),
  Expected_Heterozygosity = c(mean(sum1$Hexp), mean(sum2$Hexp), mean(sum3$Hexp), mean(sum4$Hexp), mean(sum5$Hexp), mean(sum6$Hexp), mean(sum7$Hexp)),
  Observed_Heterozygosity = c(mean(sum1$Hobs), mean(sum2$Hobs), mean(sum3$Hobs), mean(sum4$Hobs), mean(sum5$Hobs), mean(sum6$Hobs), mean(sum7$Hobs))
)

##write.table(table_Heterozygosity, "Heterozygosity.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

## PCA

# Dealing missing values
X<-tab(genind_obj, NA.method="mean")
# Running the PCA and keeping the scores for the first three PCs and then plotting the first 50 eigenvalues
pca1 <- dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 3)
pca1

png("PCA eigenvalues.png", width = 6, height = 4.5, units="in", res=1080)
barplot(pca1$eig[1:50], main = "PCA eigenvalues", col = heat.colors(50))

dev.off()

# pca individual plot
png("PCA of black pine diversity.png", width = 7, height = 5, units = "in", res = 1080)
col <- rainbow(length(levels(pop(obj))))
s.class(pca1$li, pop(genind_obj), col=colors)
title("PCA of black pine diversity \ axes 1-2")

dev.off()

# PCA individual plot with colors demonstrating genetic similarity
png("PCA of black pine_2.png", width = 6, height = 4.5, units = "in", res = 1080)
colorplot(pca1$li, pca1$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 2")
title("PCA of black pine diversity \ axes 1-2")
abline(v=0,h=0,col="grey", lty=2)

dev.off()

## DAPC

# Looking for genetic clusters
grp <- find.clusters(genind_obj, max.n.clust=10)

# A table of clusters in populations
table(pop(genind_obj), grp$grp)

# And a plot of this table
png("clusters1.png", width = 6, height = 4.5, units = "in", res = 1080)
table.value(table(pop(genind_obj), grp$grp), col.lab=paste("cluster", 1:3),
            row.lab=paste("pop", c("Dafni", "Panteleimon", "Center", "Vela", "Psaria", "Perama", "Maries")))
dev.off()

# Running the dapc

dapc1 <- dapc(genind_obj, grp$grp)


# A scatterplot of the clusters

png("clusters scatterplot.png", width = 6, height = 4.5, units = "in", res = 1080)
scatter(dapc1,1,1, col=my_colors, bg="white",
        scree.da=FALSE, legend=TRUE, solid=.4)

dev.off()

# A better one!

png("scatterplot.png", width = 6, height = 4.5, units = "in", res = 1080)
myCol <- c("darkblue","purple","green","orange","red","blue","yellow")
scatter(dapc1, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=myCol, solid=.4, cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:3))

dev.off()
# Call the dapc dunction with the genlight object and its populations
dapc.x <- dapc(genind_obj, n.pca=30,n.da=2)

# Plot the populations with the dapc function
png("dapc1.png", width = 8, height = 6, units = "in", res=1080)
scatter(dapc.x, scree.da=FALSE, cell=1.5, cex=2, bg="white", cstar=0)
title("DAPC of Black Pine Diversity Across Seven Populations")
dev.off()


## sPCA


# Using poppr - gac

library(poppr)

gac1 <- genotype_curve(obj_par, sample = 1000, quiet = TRUE)
title("Genotype accumulation curve - parents")

gac2 <- genotype_curve(obj_pro, sample = 1000, quiet = TRUE)
title("Genotype accumulation curve - progeny")

# A locus table

table1 <- locus_table(obj_par)
table1
table2 <- locus_table(obj_pro)
table2

# diversity indexes

diversity1<-poppr(obj_par)
diversity1
diversity2<-poppr(obj_pro)
diversity2

# using pegas - HW test

library(pegas)

hwe1 <- hw.test(obj_par, B = 1000)
hwe1
hwe2 <- hw.test(obj_pro, B = 1000)
hwe2

# using lattice - HW test
library("lattice")

alpha  <- 0.05
hwe21 <- hwe1
hwe21[hwe1 > alpha] <- 1
hwe22 <- hwe1
hwe22[hwe2 > alpha] <- 1


levelplot(t(hwe21))
levelplot(t(hwe22))

library("vegan")
par(mfrow=c(2,1))
tab1 <- mlg.table(obj_par, plot = FALSE)
min_sample1 <- min(rowSums(tab1))
rarecurve(tab1, sample = min_sample1, xlab = "Sample Size", ylab = "Expected MLGs")

tab2 <- mlg.table(obj_pro, plot = FALSE)
min_sample2 <- min(rowSums(tab2))
rarecurve(tab2, sample = min_sample2, xlab = "Sample Size", ylab = "Expected MLGs")


N1      <- diversity1$N      # number of samples
lambda1 <- diversity1$lambda # Simpson's index
CSim1<-(N1/(N1 - 1)) * lambda1              # Corrected Simpson's index

N2      <- diversity2$N      # number of samples
lambda2 <- diversity2$lambda # Simpson's index
CSim2<-(N2/(N2 - 1)) * lambda2              # Corrected Simpson's index

mlg.table(obj_par)
mlg.table(obj_pro)


## sPCA

# For all trees and populations

mySpca <- spca(obj, ask=FALSE, type=5, d1=0.04, d2=0.1, scannf=FALSE)

png("spatialauto2.png", width = 1000, height = 1000, res=300)
barplot(mySpca$eig,main="Eigenvalues of sPCA", col=rep(c("red","grey"),c(1,100)))
dev.off()

mySpca$as


# First deal with NAs
huhu <- tab(obj, NA.method="mean")

myGtest <- global.rtest(huhu,mySpca$lw,nperm=9999)
plot(myGtest)

myLtest <- local.rtest(huhu,mySpca$lw,nperm=9999)
plot(myLtest)

plot(mySpca)

png("colorplot2.png", width = 3000, height = 2000, res=300)
colorplot(mySpca,cex=3,main="sPCA colorplot, first global score")
dev.off()

