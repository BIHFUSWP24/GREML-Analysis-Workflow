import anndata
from sklearn.decomposition import PCA

input_file = snakemake.input["input_file"]
output_file = snakemake.output["output_file"]
dimensions = int(snakemake.wildcards["dimensions"])

# input_file = "/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/100k_snp_chr1.h5ad"
# output_file = "/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/100k_snp_chr1/PCA/8211dims.h5ad"
# dimensions = 8211

adata = anndata.read_h5ad(input_file, backed='r')
print("Data structure:")
print(adata)
print(f"Input first 5 rows and columns:")
print(adata.X[1:5, 1:5])

print(f"PCA with {dimensions} Dimensions on {adata.shape[0]} samples and {adata.shape[1]} features.")
pca = PCA(n_components=dimensions)
pca_results = pca.fit_transform(adata.X)

explained_variance = [round(var, 4) for var in pca.explained_variance_ratio_]
# print(f"Explained variance per Dimension: {explained_variance}")
# print(f"Explained variance in total: {round(sum(explained_variance),4)}")
print(f"output has shape {pca_results.shape} with first 5 rows:")
print(pca_results[1:5, :])

pca_adata = anndata.AnnData(X=pca_results, obs=adata.obs)
pca_adata.write_h5ad(output_file)
