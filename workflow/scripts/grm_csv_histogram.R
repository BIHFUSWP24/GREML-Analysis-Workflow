grm_input=snakemake@input[["grm"]]
output_file=snakemake@output[["file"]]
profile=snakemake@wildcards[["profile"]]

grm_table <- read.table(grm_input, header = FALSE, sep = ",")
grm_matrix <- as.matrix(grm_table, header = FALSE)

sorted_genetic_relationships <- sort(grm_matrix)

png(output_file, width = 1200, height = 1200, res = 150)
hist(sorted_genetic_relationships,
    breaks = "Freedman-Diaconis",
    col = "black",
    xlab = "genetic relationship",
    ylab = "number of pairs",
    main = paste("Genetic relationships Histogram of", profile),
    )
dev.off()
