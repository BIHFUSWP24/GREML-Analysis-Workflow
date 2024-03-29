grm_input=snakemake@input[["grm"]]
output_file=snakemake@output[["file"]]

grm_table <- read.csv(grm_input, header = FALSE, sep = "\t")
genetic_relationships <- grm_table[,4]
n <- grm_table[nrow(grm_table), 1]

fullGRM <- matrix(0, nrow = n, ncol = n, dimnames = list(NULL, NULL))
fullGRM[upper.tri(fullGRM, diag = TRUE)] <- genetic_relationships
fullGRM <- fullGRM + t(fullGRM) - diag(diag(fullGRM))

write.table(fullGRM, output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
