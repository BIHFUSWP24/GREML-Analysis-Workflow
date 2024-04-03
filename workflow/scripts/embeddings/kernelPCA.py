import anndata
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap
from sklearn.decomposition import KernelPCA
import matplotlib.pyplot as plt
import pandas as pd
from scipy.spatial.distance import pdist, squareform
import numpy as np
import seaborn as sns

data_path = '/sc-projects/sc-proj-dh-ag-eils-ml/genotype_data/900k_snp_chr1.h5ad'
build_path = '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/'
result_path = '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/results/'

dimensions = 2
adata = anndata.read_h5ad(data_path, backed='r')

kernel_pca = KernelPCA(n_components=dimensions, kernel='rbf')
kernel_pca_results = kernel_pca.fit_transform(adata.X)
kernel_pca_df = pd.DataFrame(kernel_pca_results, columns=['KernelPCA1', 'KernelPCA2'])
kernel_pca_df.to_csv('results/kernel_pca_results.csv', index=False)

plt.scatter(kernel_pca_results[:, 0], kernel_pca_results[:, 1])
plt.title('Kernel PCA Ergebnis')
plt.xlabel('Kernel Hauptkomponente 1')
plt.ylabel('Kernel Hauptkomponente 2')
plt.savefig(result_path + 'pca_results.png')