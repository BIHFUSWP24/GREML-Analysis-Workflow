input_file <- '/sc-projects/sc-proj-dh-ag-eils-ml/phenotype_data/continuous_phenotypes.csv'
output_file <- '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/phenotypes/phenotype_summary.tsv'

phenotype_data <- read.csv(input_file, header = TRUE, sep = ",")
ord <- order(colSums(is.na(phenotype_data)))[-1]

summary_df <- data.frame(Phenotype = character(),
    Min = numeric(),
    First_Quartile = numeric(),
    Median = numeric(),
    Mean = numeric(),
    Third_Quartile = numeric(),
    Max = numeric(),
    NA_count = integer())


for (variable in names(phenotype_data)[ord]) {
    summary_stats <- summary(phenotype_data[[variable]])
    na_count <- sum(is.na(phenotype_data[[variable]]))
    
    temp_df <- data.frame(Phenotype = variable,
        Min = summary_stats[["Min."]],
        First_Quartile = summary_stats[["1st Qu."]],
        Median = summary_stats[["Median"]],
        Mean = summary_stats[["Mean"]],
        Third_Quartile = summary_stats[["3rd Qu."]],
        Max = summary_stats[["Max."]],
        NA_count = na_count)
    summary_df <- rbind(summary_df, temp_df)
}

write.table(summary_df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
