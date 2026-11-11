#!/usr/bin/env Rscript
# Small real-data check of the shared models and dataset-specific preparation.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args[1] %in% c('kang', 'schmidt')) stop('Usage: Rscript scripts/smoke-perturb.R kang|schmidt')
suppressPackageStartupMessages({library(tidyverse); library(SingleCellExperiment); library(data.table); library(Matrix)})
options(degpaper.nebula.ncore = 1L)
source('src/common/methods.R'); source('src/common/helpers.R')
if (args[1] == 'kang') {
  source('src/kang/analysis.R')
  obj <- prepare_kang(readRDS('data/kang/Kang18_8vs8.rds'))
  cd <- as.data.frame(SummarizedExperiment::colData(obj))
  ct <- names(sort(table(cd$cluster_id), decreasing = TRUE))[1]
  ix <- which(cd$cluster_id == ct)
  group <- interaction(cd$patient_id[ix], cd$group_id[ix], drop = TRUE)
  ix <- unlist(lapply(split(ix, group), head, 20L), use.names = FALSE)
  obj <- obj[, ix]
  genes <- order(Matrix::rowMeans(SummarizedExperiment::assay(obj, 'counts')), decreasing = TRUE)[1:100]
  result <- analyze_kang(obj[genes, ])
  stopifnot(all(model_labels %in% result$method), nrow(result) > 0)
} else {
  source('src/schmidt/analysis.R')
  x <- read_schmidt('data/schmidt/Perturb-seq/data/HuTcellsCRISPRaPerturbSeq_Re-stimulated.h5Seurat')
  target <- names(sort(table(x$metadata$gene[x$metadata$gene != 'NO-TARGET']), decreasing = TRUE))[1]
  ix <- which(x$metadata$gene %in% c('NO-TARGET', target))
  group <- interaction(x$metadata$gene[ix], x$metadata$donor[ix], drop = TRUE)
  ix <- unlist(lapply(split(ix, group), head, 100L), use.names = FALSE)
  meta <- x$metadata[ix, , drop = FALSE]; mat <- x$counts[, ix, drop = FALSE]
  genes <- rownames(mat)[order(Matrix::rowMeans(mat), decreasing = TRUE)[1:100]]
  result <- analyze_schmidt_target(mat, meta, target, genes)
  stopifnot(length(result) == length(model_labels), all(vapply(result, nrow, integer(1)) > 0))
}
cat(args[1], 'smoke check completed for all six methods\n')
