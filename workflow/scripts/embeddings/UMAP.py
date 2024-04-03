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

reducer = umap.UMAP()
umap_results = reducer.fit_transform(adata.X)
umap_df = pd.DataFrame(umap_results, columns=['UMAP1', 'UMAP2'])
umap_df.to_csv('results/umap_results.csv', index=False)

plt.scatter(umap_results[:, 0], umap_results[:, 1])
plt.title('UMAP Ergebnis')
plt.xlabel('Dimension 1')
plt.ylabel('Dimension 2')
plt.savefig(result_path + 'pca_results.png')