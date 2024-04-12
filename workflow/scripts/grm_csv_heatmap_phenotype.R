suppressPackageStartupMessages(library(ComplexHeatmap))
library(colorspace)
suppressPackageStartupMessages(library(circlize))

suppressPackageStartupMessages(library(viridis))
library(hrbrthemes)
library(ggplot2)
library(reshape2)

grm_input = '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/100k_snp_chr1_pca_10dims.h5ad'
phenotype_input = '/sc-projects/sc-proj-dh-ag-eils-ml/phenotype_data/continuous_phenotypes.csv'
output_file = '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/results/phenotypes.png'
distance_method = 'correlation'

print("Reading phenotype data...")
phenotype_data <- read.csv(phenotype_input, header = TRUE, sep = ",")
print(paste0("Phenotype data dimension:(", dim(phenotype_data)[1], ",", dim(phenotype_data)[2], ")"))

print("Melting phenotype data...")
phenotype_data_long <- reshape2::melt(phenotype_data, id.vars = "eid")
names(phenotype_data_long) <- c("Sample", "Phenotype", "Value")
print(paste0("Phenotype data long dimension: (", dim(phenotype_data_long)[1], ",", dim(phenotype_data_long)[2], ") with head:"))
print(head(phenotype_data_long))

print("Plotting...")
png(output_file, width = 8000, height = 1200, res = 150)
ggplot(phenotype_data_long, aes(x=Phenotype, y=Value, fill=Phenotype)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
    theme_ipsum() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle("Boxplot") +
    xlab("")
dev.off()
