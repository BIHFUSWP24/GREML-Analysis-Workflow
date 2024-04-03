import anndata
from sklearn.decomposition import PCA
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# from sklearn.manifold import TSNE
# import umap
# from sklearn.decomposition import KernelPCA
# from scipy.spatial.distance import pdist, squareform
# import seaborn as sns

data_path = '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/100k_snp_chr1.h5ad'
csv_path = '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/pca_variances.csv'
png_path = '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/results/pca_variances.png'
adata = anndata.read_h5ad(data_path, backed='r')

pca = PCA()
pca_results = pca.fit_transform(adata.X)
variance_explained = pca.explained_variance_ratio_
variance_explained_sum = np.cumsum(variance_explained)
pca_df = pd.DataFrame({
  'Dimension': range(1, len(variance_explained)+1),
  'Variance Explained': variance_explained,
  'Cumulative Variance Explained': variance_explained_sum})
pca_df.to_csv(csv_path, index=False)

pct99 = np.argmax(variance_explained_sum >= 0.99) + 1
pct95 = np.argmax(variance_explained_sum >= 0.95) + 1

plt.plot(pca_df['Dimension'], pca_df['Cumulative Variance Explained'])
plt.xlabel('Number of Dimensions')
plt.ylabel('Cumulative Variance Explained')
plt.title('Variance Explained by PCA')
plt.axhline(y=0.99, color='r', linestyle='--')
plt.axhline(y=0.95, color='black', linestyle='--')
plt.legend(['Cumulative Variance Explained', f'99% Variance Explained ({pct99} Dimensions)', f'95% Variance Explained ({pct95} Dimensions)'], loc='lower right')
plt.savefig(png_path, dpi=300)