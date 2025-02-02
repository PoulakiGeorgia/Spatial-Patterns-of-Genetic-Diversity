# Machine-Learning

This project focuses on clustering unlabelled population genetics data. The dataset is transformated into a binary format applying one hot encoding, while also incorporating infromation on missing values and heterozygosity for each locus. Instead of clustering in high-dimensional data, an autoencoder is used to reduce the dimensionality and remove the noise. The clustering is then performed in a lower-dimensional latent space, extracted from the autoencoder, leading to a more meaningful representation of the data.  

The project is divided into three parts:
  - **One-Hot Encoding**: Data transformation into a binary format, including heterozygosity and missing values.
  - **Autoencoder**: Dimensionality reduction and noise removing.
  - **Clustering**: K-means & Hierarchical clustering implementation.
