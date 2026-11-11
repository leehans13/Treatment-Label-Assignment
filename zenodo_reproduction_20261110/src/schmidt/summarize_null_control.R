#!/usr/bin/env Rscript
# Summarize previously saved model fits; no model fitting or input data loading.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) stop('Usage: Rscript src/schmidt/summarize_null_control.R [results directory] [output directory]')
input <- if (length(args)) args[1] else 'data/schmidt/null_control'
output <- if (length(args) >= 2L) args[2] else input
source('src/common/helpers.R')

summarize_null <- function(results, test_genes, output) {
  # Each row is one gene in one randomized split. Report failures separately;
  # never turn missing p-values into non-rejections.
  rows <- lapply(seq_along(results), function(i) {
    out <- data.frame(split = i, gene = test_genes)
    for (method in model_labels) {
      fit <- results[[i]][[method]]
      out[[method]] <- if (is.null(fit)) NA_real_ else fit[match(test_genes, rownames(fit)), 'Pr(>|t|)']
    }
    out
  })
  df <- do.call(rbind, rows)
  mat <- as.matrix(df[, model_labels, drop = FALSE])
  mat[!is.finite(mat) | mat < 0 | mat > 1] <- NA_real_
  df[, model_labels] <- mat
  write.csv(df, file.path(output, 'schmidt.null.pvalues.csv'), row.names = FALSE)
  failures <- data.frame(method = model_labels, missing_or_invalid = colSums(is.na(mat)), total = nrow(mat))
  write.csv(failures, file.path(output, 'schmidt.null.failures.csv'), row.names = FALSE)
  common <- mat[complete.cases(mat), summary_order, drop = FALSE]
  if (nrow(common) < 2L) stop('Fewer than two gene/split rows have valid p-values for all methods')
  for (cutoff_name in c('bonferroni', 'nominal')) {
    cutoff <- if (cutoff_name == 'bonferroni') 0.05 / 4807 else 0.05
    sig <- common < cutoff
    overlap <- 100 * crossprod(sig) / nrow(sig)
    write.csv(round(overlap, 2), file.path(output, paste0('schmidt.null.', cutoff_name, '.pow.csv')))
  }
  sp <- suppressWarnings(cor(common, method = 'spearman'))
  write.csv(sp, file.path(output, 'schmidt.null.sp.csv'))
  # Replicate-level summaries preserve Monte Carlo variability and permit FWER calculation.
  split_summary <- do.call(rbind, lapply(seq_along(results), function(i) {
    z <- mat[df$split == i, , drop = FALSE]
    do.call(rbind, lapply(model_labels, function(method) {
      p <- z[, method]; valid <- is.finite(p)
      data.frame(split = i, method = method, tested = sum(valid), missing = sum(!valid),
        nominal_rate = if (any(valid)) mean(p[valid] < 0.05) else NA_real_,
        bonferroni_rate = if (any(valid)) mean(p[valid] < 0.05 / 4807) else NA_real_,
        any_bonferroni = if (all(valid)) any(p < 0.05 / 4807) else NA)
    }))
  }))
  write.csv(split_summary, file.path(output, 'schmidt.null.per_split.csv'), row.names = FALSE)
  message('Complete gene/split rows: ', nrow(common), ' / ', nrow(mat))
  invisible(list(pvalues = df, failures = failures, per_split = split_summary))
}

results <- readRDS(file.path(input, 'schmidt.null.results.rds'))
manifest <- readRDS(file.path(input, 'checkpoints', 'manifest.rds'))
test_genes <- manifest$extra$test_genes
if (!is.list(results) || !length(results) || !is.character(test_genes) || !length(test_genes))
  stop('Invalid saved model results or test gene list')
dir.create(output, recursive = TRUE, showWarnings = FALSE)
summarize_null(results, test_genes, output)
