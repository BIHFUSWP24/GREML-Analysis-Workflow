import anndata
from sklearn.manifold import TSNE

input_file = snakemake.input["input_file"]
output_file = snakemake.output["output_file"]
dimensions = int(snakemake.wildcards["dimensions"])
perplexity = int(snakemake.wildcards["perplexity"])

adata = anndata.read_h5ad(input_file, backed='r')
print(f"t-SNE with {dimensions} Dimensions on {adata.shape[0]} samples and {adata.shape[1]} features.")

tsne = TSNE(n_components=dimensions, perplexity=perplexity, method='exact', n_jobs=-1)
tsne_results = tsne.fit_transform(adata.X)

print(f"output has shape {tsne_results.shape} with first 5 rows:")
print(tsne_results[1:5, :])

tsne_adata = anndata.AnnData(X=tsne_results)
tsne_adata.write_h5ad(output_file)
