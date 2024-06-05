suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

# heritability_input <- "/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/results/ukb_imp_v3/analysis1_heritability_summary.tsv"
# output_path <- "/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/results/ukb_imp_v3/analysis1_summary.pdf"

heritability_input <- snakemake@input[['heritability_input']]
output_path <- snakemake@output[['summary_output']]

heritability_data <- read.csv(heritability_input, header = TRUE, check.names = FALSE, sep = "\t")
profile_names <- unique(heritability_data$Profile)
heritability_data$Profile <- factor(heritability_data$Profile, levels = profile_names)

p <- ggplot(heritability_data, aes(x = Profile, y = Heritability, fill = Category)) +
    geom_bar(stat = "identity", position = "identity") +
    geom_errorbar(aes(ymin = Heritability - SE, ymax = Heritability + SE), width = 0.2) +
    facet_wrap(~Phenotype, scales = "free_y") +
    labs(x = "Profile", y = "Heritability", fill = "Category")

ggsave(output_path, plot = p, width = 20, height = 10, dpi = 300)
