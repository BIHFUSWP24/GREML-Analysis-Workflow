grm_input <- c(
    '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/100k_snp_chr1_pca_10dims_euclidean.grm.csv',
    '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/100k_snp_chr1_pca_10dims_correlation.grm.csv',
    '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/100k_snp_chr1_pca_10dims_cosine.grm.csv')
output_file_100 <- c(
    '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/100k_snp_chr1_pca_10dims_size100_euclidean.grm.csv',
    '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/100k_snp_chr1_pca_10dims_size100_correlation.grm.csv',
    '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/100k_snp_chr1_pca_10dims_size100_cosine.grm.csv')
output_file_200 <- c(
    '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/100k_snp_chr1_pca_10dims_size200_euclidean.grm.csv',
    '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/100k_snp_chr1_pca_10dims_size200_correlation.grm.csv',
    '/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/lesi11/build/100k_snp_chr1_pca_10dims_size200_cosine.grm.csv')

resize <- function(input_file, output_file_100, output_file_200) {
    grm_table <- read.table(input_file, header = FALSE, sep = ",")
    grm_matrix <- as.matrix(grm_table, header = FALSE)
    grm_matrix_100 <- grm_matrix[1:100, 1:100]
    grm_matrix_200 <- grm_matrix[1:200, 1:200]

    write.table(grm_matrix_100, output_file_100, quote = FALSE, sep=',', col.names=FALSE, row.names = FALSE)
    write.table(grm_matrix_200, output_file_200, quote = FALSE, sep=',', col.names=FALSE, row.names = FALSE)
}

for (i in 1:length(grm_input)) {
    resize(grm_input[i], output_file_100[i], output_file_200[i])
}