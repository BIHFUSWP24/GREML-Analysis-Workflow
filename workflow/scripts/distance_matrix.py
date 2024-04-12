import anndata
import pandas as pd
from sklearn.metrics import pairwise_distances
import os

# euclidean
# cosine
# correlation
# mahalanobis


input_path = snakemake.input['input_file']
csv_path = snakemake.output['csv_output']
grm_path = snakemake.output['grm_output']
method = snakemake.wildcards['distance_method']

print("Loading data...")
adata = anndata.read_h5ad(input_path, backed='r')

df = adata.to_df()
print(f"Data has {adata.shape[0]} samples and {adata.shape[1]} features with first 5 rows and columns:")
print(df.iloc[:5, :5])

print(f"Calculating distance matrix...")
distance_matrix = pairwise_distances(df, metric=method, n_jobs=-1)
print(f"Creating distance matrix table with the {method} method...")
distance_matrix_df = pd.DataFrame(distance_matrix)
print(f"Distance matrix has shape {distance_matrix.shape} with first 5 rows and columns:")
print(distance_matrix_df.iloc[:5, :5])

print("Creating flat distance matrix table...")
n = distance_matrix.shape[0]
distance_matrix_list = [(i+1, j+1, distance_matrix[i, j]) for i in range(n) for j in range(n)]
distance_matrix_df_flat = pd.DataFrame(distance_matrix_list, columns=['patient1', 'patient2', 'distance'])
print(f"Distance matrix table has length {len(distance_matrix_df_flat['distance'])} with first 5 elements:")
print(distance_matrix_df_flat.head(5))

print("Saving distance matrix csv...")
distance_matrix_df.to_csv(csv_path, index=False, header=False)
print("Saving distance matrix grm...")
distance_matrix_df_flat.to_csv(grm_path, index=False, header=False, sep='\t')
print("Compressing distance matrix table...")
os.system(f"pigz -k -p8 {grm_path}")
