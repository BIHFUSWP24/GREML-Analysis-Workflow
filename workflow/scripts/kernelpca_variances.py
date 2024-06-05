import anndata
from sklearn.decomposition import KernelPCA
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Datenpfade
# data_input = '/sc-projects/sc-proj-dh-ag-eils-ml/genotype_data/100k_snp_all_chr.h5ad'
# csv_output = '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/kernelPCA_variances.csv'
# png_output = '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/results/kernelPCA_variances.png'

data_input = snakemake.input["data_input"] # type: ignore
csv_output = snakemake.output["csv_output"] # type: ignore
png_output = snakemake.output["png_output"] # type: ignore

# Daten laden
adata = anndata.read_h5ad(data_input, backed='r')

kpca = KernelPCA(kernel='rbf', fit_inverse_transform=True, n_components=None)  # n_components gibt die Anzahl der gewünschten Dimensionen an
kpca_results = kpca.fit_transform(adata.X)

lambdas = np.var(kpca_results, axis=0)
variance_explained = lambdas / np.sum(lambdas)
variance_explained_sum = np.cumsum(variance_explained)

pca_df = pd.DataFrame({
    'Dimension': range(1, len(variance_explained)+1),
    'Variance Explained': variance_explained,
    'Cumulative Variance Explained': variance_explained_sum
})
pca_df.to_csv(csv_output, index=False)

pct99 = np.argmax(variance_explained_sum >= 0.99) + 1
pct90 = np.argmax(variance_explained_sum >= 0.90) + 1
pct100 = np.argmax(variance_explained_sum >= 1.0) + 1

plt.plot(pca_df['Dimension'], pca_df['Cumulative Variance Explained'])
plt.xlabel('Number of Dimensions')
plt.ylabel('Cumulative Variance Explained')
plt.title('Variance Explained by KernelPCA')
plt.axhline(y=1.0, color='g', linestyle='--')
plt.axhline(y=0.99, color='r', linestyle='--')
plt.axhline(y=0.95, color='black', linestyle='--')
plt.legend(['Cumulative Variance Explained', f'100% Variance Explained ({pct100} Dimensions)' ,f'99% Variance Explained ({pct99} Dimensions)', f'90% Variance Explained ({pct90} Dimensions)'], loc='lower right')
plt.savefig(png_output, dpi=300)
