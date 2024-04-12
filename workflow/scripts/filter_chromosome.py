import anndata

data_path = '/sc-projects/sc-proj-dh-ag-eils-ml/genotype_data/100k_snp_all_chr.h5ad'
output_path = '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/100k_snp_chr1.h5ad'
output_path_tsv = '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/100k_snp_chr1.tsv'

adata = anndata.read_h5ad(data_path, backed='r')
adata_subset = adata[:, adata.var['GeneticChromosome'] == 1]
adata_subset.write_h5ad(output_path)
