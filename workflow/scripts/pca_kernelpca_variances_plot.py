import anndata
from sklearn.decomposition import PCA
from sklearn.decomposition import KernelPCA
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data_path = '/sc-projects/sc-proj-dh-ag-eils-ml/genotype_data/100k_snp_all_chr.h5ad'
png_path = '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/results/pca_kernelpca_variances.png'
adata = anndata.read_h5ad(data_path, backed='r')

# PCA
pca = PCA(n_components=None)
pca_results = pca.fit_transform(adata.X)
pca_variance_explained = pca.explained_variance_ratio_

# KernelPCA
kpca = KernelPCA(kernel='rbf', fit_inverse_transform=True, n_components=None)  # n_components gibt die Anzahl der gewünschten Dimensionen an
kpca_results = kpca.fit_transform(adata.X)

lambdas = np.var(kpca_results, axis=0)
kernelpca_variance_explained = lambdas / np.sum(lambdas)

max_len = max(len(pca_variance_explained), len(kernelpca_variance_explained))
pca_variance_explained = np.pad(pca_variance_explained, (0, max_len - len(pca_variance_explained)), 'constant')
kernelpca_variance_explained = np.pad(kernelpca_variance_explained, (0, max_len - len(kernelpca_variance_explained)), 'constant')

pca_variance_explained_sum = np.cumsum(pca_variance_explained)
kernelpca_variance_explained_sum = np.cumsum(kernelpca_variance_explained)

pca_kernelpca_df = pd.DataFrame({
    'Dimension': range(1, max_len+1),
    'PCA Variance Explained': pca_variance_explained,
    'PCA Cumulative Variance Explained': pca_variance_explained_sum,
    'KernelPCA Variance Explained': kernelpca_variance_explained,
    'KernelPCA Cumulative Variance Explained': kernelpca_variance_explained_sum
})

plt.plot(pca_kernelpca_df['Dimension'], pca_kernelpca_df['PCA Cumulative Variance Explained'])
plt.plot(pca_kernelpca_df['Dimension'], pca_kernelpca_df['KernelPCA Cumulative Variance Explained'])
plt.xlabel('Number of Dimensions')
plt.ylabel('Cumulative Variance Explained')
plt.title('Variance Explained by PCA and KernelPCA')
plt.legend(['PCA Cumulative Variance Explained', 'KernelPCA Cumulative Variance Explained'])
plt.savefig(png_path, dpi=300)