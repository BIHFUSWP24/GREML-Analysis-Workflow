import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

input_file = "/sc-projects/sc-proj-dh-ag-eils-ml/phenotype_data/continuous_phenotypes.csv"
output_file = "/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/results/phenotype_boxplot.png"

# input file format:
# Phenotype	Min	First_Quartile	Median	Mean	Third_Quartile	Max	Tenth_Percentile	Ninetieth_Percentile	NA_count
# questionaire_townsend_deprivation_index_at_recruitment_f189_0_0	-1.6040181	-0.7580106	-0.27187717	0.000127092472507085	0.59584945	3.9724765	-1.05972156	1.5205722	616

# Daten einlesen
data = pd.read_csv(input_file, header=0)

# Erstellen von zwei Subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 10), gridspec_kw={'height_ratios': [5, 1]}, sharex=True)

# Boxplot-Werte für beide Plots
data = data.drop('eid', axis=1)
list_of_arrays = [data[col].dropna().values for col in data.columns]
list_of_names = data.columns
summary = data.describe()
min_values = summary.loc['min']
max_values = summary.loc['max']

print(list_of_names)

# Erster Plot (ohne Ausreißer)
bp1 = ax1.boxplot(list_of_arrays, patch_artist=True, showfliers=False)
ax1.set_title('Boxplot of Statistical Measures for Various Phenotypes', fontsize=16)
ax1.set_ylabel('Values', fontsize=14)

# Zweiter Plot (nur Ausreißer)
ax2.boxplot(list_of_arrays,
                    boxprops={'linewidth': 0},
                    whiskerprops={'linewidth': 0},
                    capprops={'linewidth': 0},
                    flierprops={'marker': 'x', 'markersize': 4, 'alpha': 0.1},
                    medianprops={'linewidth': 0},
                    showmeans=False)
ax2.set_ylabel('Values of outliers', fontsize=14)
ax2.scatter(range(1, len(list_of_names) + 1), min_values, marker='^', color='red', zorder=3, s=16)
ax2.scatter(range(1, len(list_of_names) + 1), max_values, marker='v', color='red', zorder=3, s=16)


# X-Achse für beide Plots
for ax in [ax1, ax2]:
    ax.set_xticks(range(1, len(list_of_names) + 1))
    ax.set_xticklabels(list_of_names, rotation=45, ha='right', fontsize=8)
    ax.grid(True)

plt.tight_layout()
plt.savefig(output_file, dpi=300)
