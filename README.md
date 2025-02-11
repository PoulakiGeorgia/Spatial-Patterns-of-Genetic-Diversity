# Spatial Patterns of Genetic Diversity

## Adegenet
Genetic metrics are calculated for the whole dataset. Expected and observed heterozygosity (_He_, _Ho_), allelic richness and F-statistics are estimated. Additionally, nei's distance is computed constructing an UPGMA dendrogram. 

adegenet: a R package for the multivariate analysis of the molecular markers. 
  - **PCA**
  - **DAPC**
  - **sPCA**

## Machine Learning
It focuses on clustering unlabelled population genetics data, specifically microsatellite data with biallelic loci. The dataset is transformated into a binary format applying one hot encoding, while also incorporating infromation on missing values and heterozygosity for each locus. Instead of clustering in high-dimensional data, an autoencoder is used to reduce the dimensionality and remove the noise. The clustering is then performed in a lower-dimensional latent space, extracted from the autoencoder, leading to a more meaningful representation of the data.  

It is divided into three parts:
  - **One-Hot Encoding**: Data transformation into a binary format, including heterozygosity and missing values.
  - **Autoencoder**: Dimensionality reduction and noise removing.
  - **Clustering**: K-means & Hierarchical clustering implementation.
