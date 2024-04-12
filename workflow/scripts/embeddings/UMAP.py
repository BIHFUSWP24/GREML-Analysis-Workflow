import anndata
import umap

input_file = snakemake.input["input_file"]
output_file = snakemake.output["output_file"]
dimensions = int(snakemake.wildcards["dimensions"])
neighbors = int(snakemake.wildcards["neighbors"])

adata = anndata.read_h5ad(input_file, backed='r')

reducer = umap.UMAP(n_components=dimensions, n_neighbors=neighbors, n_jobs=-1)
umap_results = reducer.fit_transform(adata.X)

pca_adata = anndata.AnnData(X=umap_results, obs=adata.obs)
pca_adata.write_h5ad(output_file)
