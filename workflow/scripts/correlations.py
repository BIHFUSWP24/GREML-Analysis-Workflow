import anndata
import pandas as pd
import numpy as np

input_file = '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/100k_snp_chr1/kernelPCA/868dims.h5ad'
output_file = '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/100k_snp_chr1/kernelPCA/868dims/pearson_corr.grm'

adata = anndata.read_h5ad(input_file, backed='r')

corr_matrix = np.corrcoef(adata.X)

n_samples = corr_matrix.shape[0]
result = []

for i in range(n_samples):
    for j in range(i + 1):
        result.append([i + 1, j + 1, corr_matrix[i, j]])

df = pd.DataFrame(result, columns=['Sample1', 'Sample2', 'Correlation'])
print(f"Head: {df.head()}")
df.to_csv(output_file, sep='\t', index=False, header=False)
