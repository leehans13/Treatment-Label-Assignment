#!/usr/bin/env Rscript
# Count all measured genes by mean raw counts per annotated Kang cell.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) stop('Usage: Rscript src/kang/summarize_expression.R [input.rds] [output_directory]')
input <- if (length(args)) args[1] else 'data/kang/Kang18_8vs8.rds'
output <- if (length(args) >= 2L) args[2] else 'data/kang/expression_histogram'
suppressPackageStartupMessages(library(SingleCellExperiment))
source('src/kang/analysis.R')
x <- prepare_kang(readRDS(input))
means <- Matrix::rowMeans(SummarizedExperiment::assay(x, 'counts'))
if (any(!is.finite(means)) || any(means < 0)) stop('Invalid gene means')
# Lower endpoints are inclusive, upper endpoints exclusive; 1 belongs to the sixth bin.
labels <- c('< 10⁻⁵', '10⁻⁵–10⁻⁴', '10⁻⁴–10⁻³', '10⁻³–10⁻²', '10⁻²–10⁻¹', '10⁻¹–1', '> 1')
bin <- findInterval(means, c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1), left.open = FALSE) + 1L
bin[means == 1] <- 6L
counts <- tabulate(bin, nbins = length(labels))
df <- data.frame(bin = factor(labels, levels = labels), genes = counts, percent = 100 * counts / length(means))
gene_data <- data.frame(gene = rownames(x), mean_count = means, bin = labels[bin])
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write.csv(df, file.path(output, 'gene_expression_bins.csv'), row.names = FALSE)
write.csv(gene_data, file.path(output, 'gene_expression_means.csv'), row.names = FALSE)
metadata <- data.frame(genes = length(means), cells = ncol(x), zero = sum(means == 0), above_001 = sum(means > .01))
write.csv(metadata, file.path(output, 'gene_expression_metadata.csv'), row.names = FALSE)
message('Cells: ', ncol(x), '; genes: ', length(means), '; zero: ', sum(means == 0),
        '; mean > 0.01: ', sum(means > .01), ' (', round(100 * mean(means > .01), 2), '%)')
