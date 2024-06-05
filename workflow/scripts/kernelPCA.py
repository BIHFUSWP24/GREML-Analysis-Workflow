import anndata
from sklearn.decomposition import KernelPCA
from sklearn.preprocessing import StandardScaler
import pandas as pd
import os

input_file = snakemake.input["input_file"] # type: ignore
kernelpca_variances = snakemake.input["kernelpca_variances"] # type: ignore
pca_output = snakemake.output["pca_output"] # type: ignore
nms_output = snakemake.output["nms_output"] # type: ignore
kernel = snakemake.params["kernel"] # type: ignore
scale = bool(snakemake.params["scale"]) # type: ignore
variance = snakemake.params["variance"] # type: ignore
chromosomes = [int(chromosome) for chromosome in snakemake.params["chromosomes"]] # type: ignore
dimensions = snakemake.params["dimensions"] # type: ignore

if dimensions in ["All", "all"]:
    dimensions = None
    print("All dimensions")
elif dimensions is not None:
    dimensions = int(dimensions)
    print(f"Dimensions: {dimensions}")
elif variance is not None:
    variance = float(variance)
    variances_df = pd.read_csv(kernelpca_variances)
    dimensions = variances_df.loc[variances_df["Cumulative Variance Explained"] >= variance, "Dimension"]
    if not dimensions.empty:
        dimensions = int(dimensions.values[0])
    else:
        dimensions = int(variances_df.iloc[-1]["Dimension"])
    print(f"Variance: {variance} Dimensions: {dimensions}")
else:
    raise ValueError("No dimensions or variance explained specified.")

os.makedirs(os.path.dirname(pca_output), exist_ok=True)

adata = anndata.read_h5ad(input_file, backed='r')
adata = adata[:, adata.var["GeneticChromosome"].isin(chromosomes)]

print("Save non-missing SNPs...")
n_snps = adata.shape[1]
with open(nms_output, 'w') as f:
    f.write(str(n_snps))

if scale:
    print("Standardizing data...")
    scaler = StandardScaler()
    adata_matrix = scaler.fit_transform(adata.X)
else:
    adata_matrix = adata.X

print(f"Kernel PCA with {dimensions} Dimensions on {adata.shape[0]} samples and {adata.shape[1]} features.")
kernel_pca = KernelPCA(n_components=dimensions, kernel=kernel)
kernel_pca_results = kernel_pca.fit_transform(adata_matrix)

print(f"output has shape {kernel_pca_results.shape} with first 5 rows:")
print(kernel_pca_results[1:5, :])

print("Save PCA results...")
pca_adata = anndata.AnnData(X=kernel_pca_results, obs=adata.obs)
pca_adata.write_h5ad(pca_output)
