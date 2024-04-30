import anndata
from sklearn.decomposition import PCA
import os

input_file = snakemake.input["input_file"] # type: ignore
pca_output = snakemake.output["pca_output"] # type: ignore
nms_output = snakemake.output["nms_output"] # type: ignore
dimensions = int(snakemake.params["dimensions"]) # type: ignore
chromosomes = [int(chromosome) for chromosome in snakemake.params["chromosomes"]] # type: ignore

# input_file = "/sc-projects/sc-proj-dh-ag-eils-ml/genotype_data/100k_snp_all_chr.h5ad"
# pca_output = "/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/temp/PCA/8211dims.h5ad"
# non_missing_output = "/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/temp/PCA/non_missing.txt"
# dimensions = 8211
# dimensions = 10
# chromosomes = [i for i in range(1, 10)]

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

print(f"PCA with {dimensions} Dimensions on {adata.shape[0]} samples and {adata.shape[1]} features.")
pca = PCA(n_components=dimensions)
pca_results = pca.fit_transform(adata.X)

explained_variance = [round(var, 4) for var in pca.explained_variance_ratio_]
print(f"Explained variance in total: {round(sum(explained_variance),4)}")
print(f"output has shape {pca_results.shape} with first 5 rows:")
print(pca_results[1:5, :])

pca_adata = anndata.AnnData(X=pca_results, obs=adata.obs)
pca_adata.write_h5ad(pca_output)
