#!/usr/bin/env Rscript
# Run from repository root; optional arguments: number of splits, output directory.
suppressPackageStartupMessages({library(tidyverse); library(SingleCellExperiment); library(data.table); library(Matrix)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) stop('Usage: Rscript src/schmidt/run_null_control.R [splits] [output]')
n <- if (length(args)) as.integer(args[1]) else 69L
if (is.na(n) || n < 1L) stop('splits must be a positive integer')
output <- if (length(args) >= 2L) args[2] else 'data/schmidt/null_control'
input <- 'data/schmidt/Perturb-seq/data/HuTcellsCRISPRaPerturbSeq_Re-stimulated.h5Seurat'
cytokines <- 'data/schmidt/Perturb-seq/data/GO_0005125_Cytokines.txt'
options(degpaper.nebula.ncore = 1L)
source('src/common/methods.R'); source('src/common/helpers.R'); source('src/common/checkpoints.R')
source('src/schmidt/analysis.R'); source('src/schmidt/null_control.R'); source('src/schmidt/null_checkpoint_compat.R')
dir.create(output, recursive = TRUE, showWarnings = FALSE)
x <- read_schmidt(input)
genes <- schmidt_genes(x$counts, x$metadata, cytokines)
f <- hdf5r::H5File$new(input, mode = 'r')
g <- f[['meta.data/guide_id']]; guide <- g[['levels']][][g[['values']][]]; f$close_all()
selected <- x$metadata$gene == 'NO-TARGET' & x$metadata$donor %in% c('Donor1', 'Donor2')
meta <- data.frame(donor_id = x$metadata$donor[selected], guide = guide[selected], row.names = rownames(x$metadata)[selected])
meta <- meta[order(meta$donor_id), , drop = FALSE]; meta$id <- meta$donor_id
cnt <- x$counts[genes$analysis, rownames(meta), drop = FALSE]
cnt <- cnt[Matrix::rowMeans(cnt) > 0.01, , drop = FALSE]
test_genes <- unique(intersect(genes$test, rownames(cnt)))
if (anyNA(meta$guide) || length(unique(meta$donor_id)) != 2L || !length(test_genes)) stop('Invalid control metadata or test genes')
message(nrow(meta), ' control cells; ', nrow(cnt), ' analysis genes; ', length(test_genes), ' test genes')
seeds <- 20261110L + seq_len(n)
scripts <- c('src/common/methods.R', 'src/common/helpers.R', 'src/schmidt/analysis.R', 'src/schmidt/null_control.R', 'src/schmidt/run_null_control.R')
manifest <- make_manifest(c(input, cytokines), seeds, scripts, extra = list(test_genes = test_genes, analysis_genes = rownames(cnt)))
results <- checkpoint_map(seeds, output, manifest,
  compute = function(seed) run_null_models(cnt, make_null_split(meta, seed)),
  valid = function(value) is.list(value) && identical(names(value), model_labels),
  legacy = null_checkpoint_legacy)
saveRDS(results, file.path(output, 'schmidt.null.results.rds'))
message('Saved model results. Run src/schmidt/summarize_null_control.R and then src/schmidt/plot_null_control.py.')
