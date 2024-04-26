import anndata

data_path = '/sc-projects/sc-proj-dh-ag-eils-ml/genotype_data/100k_snp_all_chr.h5ad'
eid_path = '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/100k_all_eid.txt'
srid_path = '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/100k_all_srid.txt'

adata = anndata.read_h5ad(data_path, backed='r')
eids = adata.obs
print(f"eids shape: {eids.shape} with head:")
print(eids.head())

rsids = adata.var['ID']
print(f"rsids shape: {rsids.shape} with head:")
print(rsids.head())

eids.to_csv(eid_path, index=False, header=False)
rsids.to_csv(srid_path, index=False, header=False)
