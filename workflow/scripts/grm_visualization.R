library(pheatmap)
grm_input=snakemake@input[["grm"]]
output_file=snakemake@output[["file"]]
mode=snakemake@wildcards[["mode"]]

modes <- c("pheatmap", "heatmap", "histogram")
if (!(mode %in% modes)) {
  stop("mode must be one of ", paste(modes, collapse = ", "))
}

grm_table <- read.table(grm_input, header = FALSE, sep = "\t")
genetic_relationships <- grm_table[,4]
n <- grm_table[nrow(grm_table), 1]

if (mode == "pheatmap" || mode == "heatmap") {
  color_palette <- colorRampPalette(c("blue", "white", "yellow", "red"))(100)
  
  fullGRM <- matrix(0, nrow = n, ncol = n)
  fullGRM[upper.tri(fullGRM, diag = TRUE)] <- genetic_relationships
  fullGRM <- fullGRM + t(fullGRM) - diag(diag(fullGRM))
} else if (mode == "histogram") {
  sorted_genetic_relationships <- sort(genetic_relationships)
}


if (mode == "pheatmap") {
  png(output_file, width = 20000, height = 20000, res = 300)
  pheatmap(fullGRM,
      scale = "none",
      color = color_palette,
      filename = output_file,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      show_rownames = FALSE,
      show_colnames = FALSE)
  dev.off()
} else if (mode == "heatmap") {
  color_palette <- colorRampPalette(c("blue", "white", "yellow", "red"))(100)
  png(output_file, width = 20000, height = 20000, res = 300)
  heatmap(fullGRM,
          scale = "none",
          col = color_palette,
          Rowv = NA,
          Colv = NA,
          labRow = "",
          labCol = "")
  dev.off()
} else if (mode == "histogram") {
  png(output_file, width = 1200, height = 1200, res = 300)
  hist(sorted_genetic_relationships, breaks = 100, col = "blue", xlab = "genetic relationship", ylab = "number of pairs", main = "Histogram of genetic relationships", cex.main = 0.8, cex.lab = 0.8, cex.axis = 0.8, cex.sub = 0.8)
  dev.off()
}

