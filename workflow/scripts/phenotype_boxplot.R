suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

phenotype_input = '/sc-projects/sc-proj-dh-ag-eils-ml/phenotype_data/continuous_phenotypes.csv'
output_file = '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/results/phenotype_boxplot_R.png'


# Daten einlesen
print("Reading phenotype data...")
phenotype_data <- read.csv(phenotype_input, header = TRUE, sep = ",")
print(paste0("Phenotype data dimension:(", dim(phenotype_data)[1], ",", dim(phenotype_data)[2], ")"))

# Daten umformen
print("Melting phenotype data...")
phenotype_data_long <- reshape2::melt(phenotype_data, id.vars = "eid")
names(phenotype_data_long) <- c("Sample", "Phenotype", "Value")
phenotype_data_long <- phenotype_data_long[!is.na(phenotype_data_long$Value), ]
phenotype_data_long$Value <- as.numeric(as.character(phenotype_data_long$Value))
print(paste0("Phenotype data long dimension: (", dim(phenotype_data_long)[1], ",", dim(phenotype_data_long)[2], ") with head:"))
print(head(phenotype_data_long))
print(summary(phenotype_data_long$Value))

# compute the smalles and largest 10% of the data from each phenotype
percentiles_data <- phenotype_data_long %>%
    group_by(Phenotype) %>%
    summarize(
        Tenth_Percentile = quantile(Value, probs = 0.1),
        Ninetieth_Percentile = quantile(Value, probs = 0.9)
)
ylim_data <- c(min(percentiles_data$Tenth_Percentile), max(percentiles_data$Ninetieth_Percentile))

# Plotting
print("Plotting...")
png(output_file, width = 3200, height = 2000, res = 300)
ggplot(phenotype_data_long, aes(x = Phenotype, y = Value)) +
    geom_boxplot(outlier.shape = NA) +
    labs(x = "", y = "Phenotype Value", title = "Phenotype Boxplot") +
    theme_minimal() +
    scale_y_continuous(limits = ylim_data) +
    coord_cartesian(ylim = ylim_data) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1, size = 5))
dev.off()
