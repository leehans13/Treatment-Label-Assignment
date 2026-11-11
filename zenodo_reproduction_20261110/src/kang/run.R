#!/usr/bin/env Rscript
# Run from the repository root: Rscript src/kang/run.R [input.rds] [output_directory]
suppressPackageStartupMessages({
  library(tidyverse)
  library(SingleCellExperiment)
  library(data.table)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) stop("Usage: Rscript src/kang/run.R [input.rds] [output_directory]")
input <- if (length(args) >= 1L) args[1] else "data/kang/Kang18_8vs8.rds"
output <- if (length(args) >= 2L) args[2] else "data/kang/results"
source("src/common/methods.R")
source("src/common/helpers.R")
source("src/common/checkpoints.R")
source("src/kang/analysis.R")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
obj <- prepare_kang(readRDS(input))
message("Input: ", input, "; ", nrow(obj), " genes, ", ncol(obj), " cells")
units <- unique(as.character(SummarizedExperiment::colData(obj)$cluster_id))
manifest <- make_manifest(input, units, c("src/common/methods.R", "src/common/helpers.R", "src/kang/analysis.R"))
result <- analyze_kang(obj, function(units, compute) {
  checkpoint_map(units, output, manifest, compute,
                 valid = function(x) is.data.frame(x) && all(model_labels %in% x$method))
})
saveRDS(result, file.path(output, "kang.result.rds"))
summarize_kang(result, output)
message("Wrote results to ", output)
