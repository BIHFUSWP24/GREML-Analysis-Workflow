suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(colorspace))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(forcats))

set.seed(123)

grm_input=snakemake@input[["grm"]]
grmID_input=snakemake@input[["grmID"]]
phenotype_inputs=snakemake@input[["phenotype_inputs"]]
phenotype_names=snakemake@params[["phenotype_names"]]
output_file=snakemake@output[["file"]]
profile=snakemake@wildcards[["profile"]]

# grm_input <- "/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/ukb_imp_v3/grm/pca/PCA100%.grm.gz"
# grmID_input <- "/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/ukb_imp_v3/grm/pca/PCA100%.grm.id"
# phenotype_inputs <- c("/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/temp/measurements_body_mass_index_bmi_f21001_0_0.phen")
# phenotype_names <- c("BMI")
# output_file <- "/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/results/ukb_imp_v3/heatmap/PCA100%.png"
# profile <- "PCA100%"

output_file <- gsub("%", "%%", output_file)

print("Reading GRM file...")
grm_table <- read.csv(grm_input, header = FALSE, sep = "\t")
genetic_relationships <- grm_table[, 4]
n <- grm_table[nrow(grm_table), 1]

print("Creating GRM...")
fullGRM <- matrix(0, nrow = n, ncol = n, dimnames = list(NULL, NULL))
fullGRM[upper.tri(fullGRM, diag = TRUE)] <- genetic_relationships
fullGRM <- fullGRM + t(fullGRM) - diag(diag(fullGRM))

min_off_diag_val <- min(fullGRM[upper.tri(fullGRM, diag = FALSE)], na.rm = TRUE)
max_off_diag_val <- max(fullGRM[upper.tri(fullGRM, diag = FALSE)], na.rm = TRUE)
max_full_val <- max(fullGRM, na.rm = TRUE)

print(paste("GRM dimensions:", dim(fullGRM)))
print(paste("min offdiag:", min_off_diag_val, "max offdiag:", max_off_diag_val, "max:", max_full_val))

print("Creating GRM Annotation...")
print(paste("   Reading GRM ID data from", grmID_input, "..."))
grmId_data <- read.csv(grmID_input, sep = "\t", header = FALSE)
eids <- grmId_data[,1]
eids_df <- data.frame(FID = eids)

for (i in 1:length(phenotype_inputs)) {
    print(paste("   Reading phenotype data from", phenotype_inputs[i], "..."))
    phenotype_data <- read.csv(phenotype_inputs[i], sep = " ", header = FALSE)
    colnames(phenotype_data) <- c("FID", "IID", phenotype_names[i])
    phenotype_data <- phenotype_data[, c("FID", phenotype_names[i])]
    unique_values <- unique(phenotype_data[[phenotype_names[i]]])
    if (length(unique_values) == 2) {
        phenotype_data[[phenotype_names[i]]] <- factor(phenotype_data[[phenotype_names[i]]], levels = unique_values)
        if (all(unique_values %in% c(0, 1))) {
            phenotype_data[[phenotype_names[i]]] <- fct_recode(phenotype_data[[phenotype_names[i]]], `False` = "0", `True` = "1")
        }
    }
    eids_df <- merge(x = eids_df, y = phenotype_data, by = "FID", all.x = TRUE)
}

color_list <- qualitative_hcl(length(phenotype_names), palette = "Dark 2")
color_ramps <- list()
for (i in 1:length(phenotype_names)) {
    column_data <- eids_df[[phenotype_names[i]]]
    if (is.factor(column_data)) {
        levels_count <- length(levels(column_data))
        color_ramps[[phenotype_names[i]]] <- setNames(c("#FFFFFFFF", color_list[i]), levels(column_data))
    } else {
        min_val <- min(na.omit(column_data), na.rm = TRUE)
        max_val <- max(na.omit(column_data), na.rm = TRUE)
        color_ramps[[phenotype_names[i]]] <- colorRamp2(c(min_val, max_val), c("white", color_list[i]))
    }
}

print("Summary of eids_df:")
print(summary(eids_df))
print(paste("Dimensions of eids_df:", nrow(eids_df), "eids and", ncol(eids_df), "phenotypes"))
print("Head of eids_df:")
print(head(eids_df))
print("Color ramps:")
print(color_ramps)

pheno_anno <- rowAnnotation(df = eids_df[,-1],
    col = color_ramps,
    na_col = "white",
    annotation_legend_param = list(
        title_position = "leftcenter-rot",
        legend_height = unit(1.5, "in"),
        title_gp = gpar(fontsize = 8),
        labels_gp = gpar(fontsize = 8)
    )
)

print("Creating heatmap...")
if (min_off_diag_val > 0 || max_off_diag_val < 0) {
    color_ramp <- colorRamp2(c(min_off_diag_val, 0, max_off_diag_val, max_full_val), c("red", "white", "blue", "black"))
} else if (-min_off_diag_val > max_off_diag_val) {
    color_ramp <- colorRamp2(c(min_off_diag_val, -max_off_diag_val, 0, max_off_diag_val, max_full_val), c("red", "orange", "white", "blue", "black"))
} else {
    color_ramp <- colorRamp2(c(min_off_diag_val, 0, -min_off_diag_val, max_off_diag_val, max_full_val), c("red", "white", "blue", "darkmagenta", "black"))
}

heatmap <- Heatmap(fullGRM,
    use_raster = TRUE,
    raster_quality = 10,
    name = "Genetic Relationship Matrix",
    width = unit(12, "in"),
    height = unit(12, "in"),
    col = color_ramp,
    column_title = paste(profile, "GRM heatmap"),
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    show_column_dend = FALSE,
    show_row_dend = FALSE,
    border_gp = gpar(col = "black", lty = 0),
    right_annotation = pheno_anno,
    heatmap_legend_param = list(
        title = "GR",
        title_position = "leftcenter-rot",
        legend_height = unit(1.5, "in"),
        title_gp = gpar(fontsize = 8),
        labels_gp = gpar(fontsize = 8)
    ),)

print("Saving heatmap...")
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
png(output_file, width=20, height=16, units="in", res=834)
draw(heatmap)
dev.off()
