import pandas as pd

input_file = snakemake.input["input_file"]
output_file = snakemake.output["output_file"]
selected_phenotype = snakemake.wildcards["phenotype"]
mode = snakemake.wildcards["mode"]

if mode == 'original':
    phenotypes = pd.read_feather(input_file)
elif mode == 'normalized':
    phenotypes = pd.read_csv(input_file)
else:
    raise ValueError(f"Unknown mode {mode}.")

phenotypes.index = phenotypes['eid']
print(f"Dimensions of the phenotype data: {phenotypes.shape}")
if selected_phenotype not in phenotypes.columns:
    print(f"Phenotypes available: {phenotypes.columns}")
    print(f"Phenotype selected: {selected_phenotype}")
    raise ValueError(f"Phenotype {selected_phenotype} not found in {input_file}.")

print("1")
phenotype_values = phenotypes[selected_phenotype].dropna()
print("2")
results = pd.DataFrame({'Family ID': phenotype_values.index, 'Individual ID': phenotype_values.index, 'Phenotype': phenotype_values})
print("3")

results.to_csv(output_file, sep=' ', index=False, header=False)
