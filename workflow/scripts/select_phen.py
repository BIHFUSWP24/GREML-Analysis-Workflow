import pandas as pd

# output format: Family ID, Individual ID, Phenotype

input_file = '/sc-projects/sc-proj-dh-ag-eils-ml/phenotype_data/continuous_phenotypes.csv'
output_file = '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/phenotypes/hip_circumference_f49_0_0.phen'
selected_phenotype = 'measurements_hip_circumference_f49_0_0'

phenotypes = pd.read_csv(input_file, index_col=0)
print(f"Dimensions of the phenotype data: {phenotypes.shape}")
if selected_phenotype not in phenotypes.columns:
    raise ValueError(f"Phenotype {selected_phenotype} not found in {input_file}.")

phenotype_values = phenotypes[selected_phenotype].dropna()
results = pd.DataFrame({'Family ID': phenotype_values.index, 'Individual ID': phenotype_values.index, 'Phenotype': phenotype_values})

results.to_csv(output_file, sep=' ', index=False, header=False)
