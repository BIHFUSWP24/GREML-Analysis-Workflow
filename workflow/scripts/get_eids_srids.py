import anndata

input_file = snakemake.input['input_file'] # type: ignore
eid_output = snakemake.output['eid_output'] # type: ignore
srid_output = snakemake.output['srid_output'] # type: ignore

adata = anndata.read_h5ad(input_file, backed='r')
eids = adata.obs
print(f"eids shape: {eids.shape} with head:")
print(eids.head())

rsids = adata.var['ID']
print(f"rsids shape: {rsids.shape} with head:")
print(rsids.head())

eids.to_csv(eid_output, index=False, header=False)
rsids.to_csv(srid_output, index=False, header=False)
