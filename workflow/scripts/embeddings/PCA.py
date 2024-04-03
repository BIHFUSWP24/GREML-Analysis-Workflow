import anndata
from sklearn.decomposition import PCA
import pandas as pd

data_path = '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/100k_snp_chr1.h5ad'
build_path = '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/'
variances_file = '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/pca_variances.csv'

dimensions = 10
adata = anndata.read_h5ad(data_path, backed='r')
print(f"PCA with {dimensions} Dimensions on {adata.shape[0]} samples and {adata.shape[1]} features.")

pca = PCA(n_components=dimensions)
pca_results = pca.fit_transform(adata.X)

explained_variance = [round(var, 4) for var in pca.explained_variance_ratio_]
print(f"Explained variance per Dimension: {explained_variance}")
print(f"Explained variance in total: {round(sum(explained_variance),4)}")
print(f"output has shape {pca_results.shape} with first 5 rows:")
print(pca_results[1:5, :])

pca_adata = anndata.AnnData(X=pca_results)
pca_adata.write_h5ad(f'{build_path}100k_snp_chr1_pca_{dimensions}dims.h5ad')
