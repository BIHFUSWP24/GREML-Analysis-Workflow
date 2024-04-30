library(tidyr)
library(ggplot2)
library(dplyr)

# input_path <- "/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/results/ukb_imp_v3/conmpare_greml_summary.tsv"
# output_path <- "/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/results/ukb_imp_v3/conmpare_greml_summary.pdf"

input_path=snakemake@input[['summary_input']]
output_path=snakemake@output[['summary_output']]

data <- read.csv(input_path, header = TRUE, check.names = FALSE, sep = "\t")

data_long <- data %>%
    pivot_longer(cols = -Phenotype, names_to = "method", values_to = "value") %>%
    mutate(method = factor(method, levels = unique(method))) %>%
    mutate(group = ifelse(grepl("kernel", method), "kernelPCA", ifelse(grepl("PCA", method), "PCA", "none"))) %>%
    mutate(dims = factor(sub("^kernelPCA|^PCA", "", method), levels = unique(sub("^kernelPCA|^PCA", "", method))))
    
print(data_long)

# Balkendiagramm für jeden Phänotyp
p <- ggplot(data_long, aes(x = dims, y = value, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~Phenotype, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "GREML results for each Method per Phenotype",
       x = "Method",
       y = "GREML Result")

ggsave(output_path, plot = p, width = 20, height = 10, dpi = 300)
