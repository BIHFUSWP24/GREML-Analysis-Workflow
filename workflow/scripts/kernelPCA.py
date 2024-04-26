import anndata
from sklearn.decomposition import KernelPCA
from sklearn.preprocessing import StandardScaler
import os

input_file = snakemake.input["input_file"]
output_file = snakemake.output["output_file"]
dimensions = int(snakemake.params["dimensions"])
chromosomes = [int(chromosome) for chromosome in snakemake.params["chromosomes"]]

os.makedirs(os.path.dirname(output_file), exist_ok=True)

adata = anndata.read_h5ad(input_file, backed='r')
adata = adata[:, adata.var["GeneticChromosome"].isin(chromosomes)]

# is_not_nan_matrix = ~np.isnan(adata.X)
# with open(non_missing_output, "w") as f:
#     for i in range(adata.shape[0]):
#         for j in range(i+1):
#             no_nans = np.sum(np.logical_and(is_not_nan_matrix[i, :], is_not_nan_matrix[j, :]))
#             f.write(f"{i} {j} {no_nans}\n")

scaler = StandardScaler()
scaled_adata = scaler.fit_transform(adata.X)

kernel_pca = KernelPCA(n_components=dimensions, kernel='rbf')
kernel_pca_results = kernel_pca.fit_transform(scaled_adata)

pca_adata = anndata.AnnData(X=kernel_pca_results, obs=adata.obs)
pca_adata.write_h5ad(output_file)
