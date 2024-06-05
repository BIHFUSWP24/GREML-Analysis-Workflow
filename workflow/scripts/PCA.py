import anndata
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import os
import pandas as pd

input_file = snakemake.input["input_file"] # type: ignore
pca_variances = snakemake.input["pca_variances"] # type: ignore
pca_output = snakemake.output["pca_output"] # type: ignore
nms_output = snakemake.output["nms_output"] # type: ignore
chromosomes = [int(chromosome) for chromosome in snakemake.params["chromosomes"]] # type: ignore
dimensions = snakemake.params["dimensions"] # type: ignore
variance = snakemake.params["variance"] # type: ignore
scale = bool(snakemake.params["scale"]) # type: ignore

if dimensions in ["All", "all"]:
    dimensions = None
    print("All dimensions")
elif dimensions is not None:
    dimensions = int(dimensions)
    print(f"Dimensions: {dimensions}")
elif variance is not None:
    variance = float(variance)
    variances_df = pd.read_csv(pca_variances)
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

print("Data structure:")
print(adata)
print("Contains chromosomes:")
print(adata.var["GeneticChromosome"].unique())

print(f"Filtering for chromosomes {chromosomes}")
adata = adata[:, adata.var["GeneticChromosome"].isin(chromosomes)]
print(f"Filtered data structure:")
print(adata)

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

print(f"PCA with {dimensions} Dimensions on {adata.shape[0]} samples and {adata.shape[1]} features.")
pca = PCA(n_components=dimensions)
pca_results = pca.fit_transform(adata_matrix)

explained_variance = [round(var, 4) for var in pca.explained_variance_ratio_]
print(f"Explained variance in total: {round(sum(explained_variance),4)}")
print(f"output has shape {pca_results.shape} with first 5 rows:")
print(pca_results[1:5, :])

pca_adata = anndata.AnnData(X=pca_results, obs=adata.obs)
pca_adata.write_h5ad(pca_output)
