suppressPackageStartupMessages(library(ComplexHeatmap))
library(colorspace)
suppressPackageStartupMessages(library(circlize))

grm_input=snakemake@input[["grm"]]
output_file=snakemake@output[["file"]]
distance_method=snakemake@wildcards[["distance_method"]]

grm_table <- read.table(grm_input, header = FALSE, sep = ",")
grm_matrix <- as.matrix(grm_table, header = FALSE)
rownames(grm_matrix) <- rep(" ", dim(grm_matrix)[1])
colnames(grm_matrix) <- rep(" ", dim(grm_matrix)[2])

min_val <- min(grm_matrix, na.rm = TRUE)
max_val <- max(grm_matrix, na.rm = TRUE)
val = max(abs(min_val), abs(max_val))

color_ramp <- colorRamp2(c(-val, 0, val), c("red", "white", "blue"))

heatmap <- Heatmap(grm_matrix,
    raster_quality = 5,
    name = "Genetic Relationship",
    width = unit(15, "cm"),
    height = unit(15, "cm"),
    col = color_ramp,
    column_title = paste("Genetic relationships heatmap with", distance_method, "distance method"),
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    )

png(output_file, width=19, height=22, units="cm", res=1200)
draw(heatmap,
    heatmap_legend_side = "bottom",
    )
dev.off()
