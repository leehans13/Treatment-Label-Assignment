# Shared model execution and notebook-style summaries.
model_labels <- c("NB GLMM", "glmGamPoi (Pb)", "edgeR (Pb)", "limma (Pb)", "glmGamPoi (cell)", "fixest (cell)")
summary_order <- c("glmGamPoi (Pb)", "edgeR (Pb)", "limma (Pb)", "glmGamPoi (cell)", "fixest (cell)", "NB GLMM")
run_models <- function(counts, metadata, fixed_effect) {
  fun <- list(nebula.mult, glmgp.mult, edger.mult, limma.mult, glmgp.cell.mult, fixest.mult)
  result <- lapply(seq_along(fun), function(i) {
    message("  ", model_labels[i])
    form <- if (i == 6L) stats::as.formula(paste("y ~ tx_cell |", fixed_effect)) else ~ tx_cell
    mat <- if (i == 6L) Matrix::t(counts) else counts
    fun[[i]](form, mat, metadata)
  })
  names(result) <- model_labels
  result
}
write_summaries <- function(df, cutoff, output, prefix) {
  mat <- as.matrix(df[, -1, drop = FALSE])
  mat[!is.finite(mat[, 1]), ] <- 1
  sig <- mat < cutoff
  pow <- round(crossprod(sig) / nrow(sig) * 100, 2)
  sp <- stats::cor(mat, method = "spearman")
  pow <- pow[summary_order, summary_order, drop = FALSE]
  sp <- sp[summary_order, summary_order, drop = FALSE]
  utils::write.csv(pow, file.path(output, paste0(prefix, ".pow.csv")))
  utils::write.csv(sp, file.path(output, paste0(prefix, ".sp.csv")))
  invisible(list(power = pow, spearman = sp))
}
