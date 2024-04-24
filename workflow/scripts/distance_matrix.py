import anndata
import pandas as pd
from sklearn.metrics import pairwise_distances
import os
import sys

log_path = snakemake.log[0]
sys.stdout = open(log_path, 'w')
sys.stderr = sys.stdout

input_path = snakemake.input['input_file']
grm_path = snakemake.output['grm_output']
method = snakemake.wildcards['distance_method']

print("Loading data...")
adata = anndata.read_h5ad(input_path, backed='r')

adata_df = adata.to_df()
print(f"Data has {adata.shape[0]} samples and {adata.shape[1]} features with first 5 rows and columns:")
print(adata_df.iloc[:5, :5])

print(f"Calculating distance matrix...")
if method == 'pearson' or method == 'spearman':
    adata_df = adata_df.T
    distance_df = adata_df.corr(method=method)
    distance_matrix = distance_df.values
else:
    distance_matrix = pairwise_distances(adata_df, metric=method, n_jobs=-1)

print(f"Distance matrix has shape {distance_matrix.shape} with first 5 rows and columns:")
print(distance_matrix[:5, :5])

print("Creating flat distance matrix table...")
n = distance_matrix.shape[0]
distance_matrix_list = [(i+1, j+1, 0, distance_matrix[i, j]) for i in range(n) for j in range(i+1)]
distance_matrix_df_flat = pd.DataFrame(distance_matrix_list, columns=['patient1', 'patient2', 'non-missing SNPs', 'distance'])
print(f"Distance matrix table has length {len(distance_matrix_df_flat['distance'])} with first 5 elements:")
print(distance_matrix_df_flat.head(5))

print("Saving distance matrix grm...")
distance_matrix_df_flat.to_csv(grm_path, index=False, header=False, sep='\t')
print("Compressing distance matrix table...")
os.system(f"pigz -k -p8 {grm_path}")
print("Done!")

sys.stdout.close()
sys.stdout = sys.__stdout__
sys.stderr = sys.__stderr__
