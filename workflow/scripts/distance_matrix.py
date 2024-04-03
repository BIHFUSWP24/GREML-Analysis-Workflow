import anndata
import pandas as pd
from sklearn.metrics import pairwise_distances
import os

input_name = '100k_snp_chr1_pca_10dims'
method = 'euclidean'

input_path = f'/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/{input_name}.h5ad'
csv_path = f'/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/{input_name}_{method}.csv'
grm_path = f'/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/{input_name}_{method}.grm'

adata = anndata.read_h5ad(input_path, backed='r')

print(f"Calculating distance matrix for {adata.shape[0]} samples and {adata.shape[1]} features...")
df = pd.DataFrame(adata.X)
distance_matrix = pairwise_distances(df, metric=method, n_jobs=-1)
print("Creating distance matrix table...")
distance_matrix_df = pd.DataFrame(distance_matrix)
print(f"Distance matrix has shape {distance_matrix.shape} with first 5 rows and columns:")
print(distance_matrix_df.iloc[:5, :5])

print("Creating flat distance matrix table...")
n = distance_matrix.shape[0]
distance_matrix_list = [(i+1, j+1, distance_matrix[i, j]) for i in range(n) for j in range(n)]
distance_matrix_df_flat = pd.DataFrame(distance_matrix_list, columns=['patient1', 'patient2', 'distance'])
print(f"Distance matrix table has length {len(distance_matrix_df['distance'])} with first 5 elements:")
print(distance_matrix_df.head(5))

print("Saving distance matrix tables...")
distance_matrix_df.to_csv(csv_path, index=False, header=False)
distance_matrix_df_flat.to_csv(grm_path, index=False, header=False, sep='\t')
print("Compressing distance matrix table...")
os.system(f"pigz -k -p8 {grm_path}")
