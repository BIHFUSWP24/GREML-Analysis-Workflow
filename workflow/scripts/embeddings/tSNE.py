import anndata
from sklearn.manifold import TSNE

data_path = '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/100k_snp_chr1.h5ad'
build_path = '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/'

dimensions = 10
adata = anndata.read_h5ad(data_path, backed='r')
print(f"t-SNE with {dimensions} Dimensions on {adata.shape[0]} samples and {adata.shape[1]} features.")

tsne = TSNE(n_components=dimensions, method='exact', n_jobs=-1)
tsne_results = tsne.fit_transform(adata.X)

print(f"output has shape {tsne_results.shape} with first 5 rows:")
print(tsne_results[1:5, :])

tsne_adata = anndata.AnnData(X=tsne_results)
tsne_adata.write_h5ad(f'{build_path}100k_snp_chr1_tsne_{dimensions}dims.h5ad')
