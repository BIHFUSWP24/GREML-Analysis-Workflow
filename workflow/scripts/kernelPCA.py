import anndata
from sklearn.decomposition import KernelPCA
from sklearn.preprocessing import StandardScaler
import os

input_file = snakemake.input["input_file"] # type: ignore
pca_output = snakemake.output["pca_output"] # type: ignore
nms_output = snakemake.output["nms_output"] # type: ignore
dimensions = int(snakemake.params["dimensions"]) # type: ignore
chromosomes = [int(chromosome) for chromosome in snakemake.params["chromosomes"]] # type: ignore

os.makedirs(os.path.dirname(pca_output), exist_ok=True)

adata = anndata.read_h5ad(input_file, backed='r')
adata = adata[:, adata.var["GeneticChromosome"].isin(chromosomes)]

print("Save non-missing SNPs...")
n_snps = adata.shape[1]
with open(nms_output, 'w') as f:
    f.write(str(n_snps))

scaler = StandardScaler()
scaled_adata = scaler.fit_transform(adata.X)

kernel_pca = KernelPCA(n_components=dimensions, kernel='rbf')
kernel_pca_results = kernel_pca.fit_transform(scaled_adata)

pca_adata = anndata.AnnData(X=kernel_pca_results, obs=adata.obs)
pca_adata.write_h5ad(pca_output)
