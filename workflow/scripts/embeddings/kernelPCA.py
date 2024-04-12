import anndata
from sklearn.decomposition import KernelPCA
import os

input_file = snakemake.input["input_file"]
output_file = snakemake.output["output_file"]
dimensions = int(snakemake.wildcards["dims"])

adata = anndata.read_h5ad(input_file, backed='r')

kernel_pca = KernelPCA(n_components=dimensions, kernel='rbf')
kernel_pca_results = kernel_pca.fit_transform(adata.X)

os.makedirs(os.path.dirname(output_file), exist_ok=True)
pca_adata = anndata.AnnData(X=kernel_pca_results, obs=adata.obs)
pca_adata.write_h5ad(output_file)
