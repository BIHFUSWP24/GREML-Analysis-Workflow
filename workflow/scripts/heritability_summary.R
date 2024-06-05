library(tidyr)
library(dplyr)

input_files <- snakemake@input[["hsq"]]
output_file <- snakemake@output[["file"]]
methods <- snakemake@params[["methods"]]
n_phenotypes <- length(input_files)/length(methods)
method_list <- rep(methods, each = n_phenotypes)

file_names <- basename(input_files)

parts <- strsplit(file_names, ".", fixed = TRUE)
profiles <- sapply(parts, function(x) { return(x[1]) })
phenotype_list <- sapply(parts, function(x) { return(x[2]) })

heritabilities <- sapply(input_files, function(file) { 
    if (file.size(file) == 0) {
        "Error"
    } else {
        data <- read.delim(file, header = TRUE, sep = "\t", fill = TRUE)
        data[["Variance"]][4]}})

standard_errors <- sapply(input_files, function(file) { 
    if (file.size(file) == 0) {
        "Error"
    } else {
        data <- read.delim(file, header = TRUE, sep = "\t", fill = TRUE)
        data[["SE"]][4]}})

df <- data.frame(Profile=profiles, Phenotype=phenotype_list, Category=method_list, Heritability=heritabilities, SE=standard_errors, stringsAsFactors = FALSE)
rownames(df) <- NULL

write.table(df, output_file, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
