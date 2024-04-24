library(tidyr)

input_files <- snakemake@input[["hsq"]]
output_file <- snakemake@output[["file"]]
phenotypes <- snakemake@params[["phenotypes"]]
file_names <- basename(input_files)

get_categories <- function(original_string, delete_strings) {
    pattern <- paste(delete_strings, collapse="|")
    gsub(pattern, "", original_string)
}

delete <- c(paste0(phenotypes, "_"), ".hsq")
categories <- sapply(file_names, get_categories, delete)

print("Categories:")
print(categories)

heritabilities <- sapply(input_files, function(file) {
    data <- read.delim(file, header = TRUE, sep = "\t", fill = TRUE)
    data$Variance[4]
})

file_names <- gsub(".hsq", "", file_names)
for (category in categories) {
    file_names <- gsub(paste0("_", category), "", file_names)
}

df <- data.frame(Phenotype=file_names, Category=categories, Heritability=heritabilities)
row.names(df) <- NULL

df_wide <- pivot_wider(df, names_from = Category, values_from = Heritability)

write.table(df_wide, output_file, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
