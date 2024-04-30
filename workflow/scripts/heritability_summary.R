library(tidyr)
library(dplyr)

input_files <- snakemake@input[["hsq"]]
output_file <- snakemake@output[["file"]]

file_names <- basename(input_files)

print("Input files:")
print(input_files)

parts <- strsplit(file_names, ".", fixed = TRUE)
categories <- sapply(parts, function(x) { return(x[1]) })
phenotype_list <- sapply(parts, function(x) { return(x[2]) })

heritabilities <- sapply(input_files, function(file) { 
    if (file.size(file) == 0) {
        "Error"
    } else {
        data <- read.delim(file, header = TRUE, sep = "\t", fill = TRUE)
        data$Variance[4]}})

df <- data.frame(Phenotype=phenotype_list, Category=categories, Heritability=heritabilities)

df_wide <- pivot_wider(df, names_from = Category, values_from = Heritability) %>% arrange(Phenotype)
print("Wide format:")
print(df_wide)

write.table(df_wide, output_file, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
