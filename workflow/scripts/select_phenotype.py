import pandas as pd

input_file = snakemake.input["input_file"]
output_file = snakemake.output["output_file"]
selected_phenotype = snakemake.wildcards["phenotype"]

if input_file.endswith('.feather'):
    phenotypes = pd.read_feather(input_file)
elif input_file.endswith('.csv'):
    phenotypes = pd.read_csv(input_file)
else:
    raise ValueError(f"Can't open {input_file}.")

phenotypes.index = phenotypes['eid']
if selected_phenotype not in phenotypes.columns:
    print(f"Can't find phenotype {selected_phenotype} in {input_file}.")
    print(f"Phenotypes available: {phenotypes.columns}")
    raise ValueError(f"Phenotype {selected_phenotype} not found in {input_file}.")

phenotype_values = phenotypes[selected_phenotype].dropna()
results = pd.DataFrame({'Family ID': phenotype_values.index, 'Individual ID': phenotype_values.index, 'Phenotype': phenotype_values})

results.to_csv(output_file, sep=' ', index=False, header=False)
