# Randomized negative controls for the Schmidt re-stimulated dataset.
make_null_split <- function(metadata, seed) {
  set.seed(seed)
  metadata$tx_cell <- NA_integer_
  # Balance labels within donor and guide to avoid confounding either with label.
  strata <- interaction(metadata$donor_id, metadata$guide, drop = TRUE)
  for (idx in split(seq_len(nrow(metadata)), strata)) {
    n1 <- length(idx) %/% 2L
    if (length(idx) %% 2L) n1 <- n1 + sample.int(2L, 1L) - 1L
    metadata$tx_cell[idx] <- sample(c(rep.int(1L, n1), rep.int(0L, length(idx) - n1)))
  }
  metadata
}

run_null_models <- function(counts, metadata) {
  fun <- list(nebula.mult, glmgp.mult, edger.mult, limma.mult, glmgp.cell.mult, fixest.mult)
  result <- lapply(seq_along(fun), function(i) {
    message('  ', model_labels[i])
    tryCatch({
      form <- if (i == 6L) y ~ tx_cell | donor_id else ~ tx_cell
      mat <- if (i == 6L) Matrix::t(counts) else counts
      as.data.frame(fun[[i]](form, mat, metadata))
    }, error = function(e) { message('  FAILED: ', conditionMessage(e)); NULL })
  })
  names(result) <- model_labels
  result
}

