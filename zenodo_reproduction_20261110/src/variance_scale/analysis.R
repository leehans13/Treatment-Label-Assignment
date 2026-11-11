# One null draw and the same NEBULA model used in the principal benchmark.
simulate_fit <- function(scenario, n, m, seed) {
  set.seed(seed)
  id <- rep(seq_len(n), each = m)
  alpha <- rnorm(n, 0, 0.5)
  treatment <- if (scenario == 1L) rep(sample(rep(0:1, each = n / 2)), each = m) else rbinom(n * m, 1, 0.5)
  y <- rnbinom(n * m, mu = exp(alpha[id]), size = 2)
  count <- matrix(y, nrow = 1, dimnames = list('gene', NULL))
  pred <- model.matrix(~ treatment)
  warnings <- character()
  result <- tryCatch(withCallingHandlers(
    nebula::nebula(count, id, pred = pred, offset = rep(1, n * m), cpc = 0, mincp = 0, ncore = 1),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning') }),
    error = function(e) e)
  record <- data.frame(scenario = scenario, n = n, m = m, seed = seed,
    estimate = NA_real_, se = NA_real_, p = NA_real_, subject_variance = NA_real_, cell_dispersion = NA_real_,
    convergence = NA_real_, algorithm = NA_character_, boundary = NA, valid = FALSE,
    diagnostic = paste(unique(warnings), collapse = ' | '))
  if (inherits(result, 'error')) { record$diagnostic <- conditionMessage(result); return(record) }
  s <- result$summary
  record$estimate <- s$logFC_treatment[1]; record$se <- s$se_treatment[1]
  record$p <- s$p_treatment[1]; record$convergence <- result$convergence[1]
  record$algorithm <- result$algorithm[1]
  record$subject_variance <- result$overdispersion$Subject[1]
  record$cell_dispersion <- result$overdispersion$Cell[1]
  # NEBULA defaults bound subject overdispersion below by 1e-4.
  record$boundary <- if (is.finite(record$subject_variance)) record$subject_variance <= 1.0001e-4 else NA
  record$valid <- is.finite(record$convergence) && record$convergence > -20 &&
    is.finite(record$estimate) && is.finite(record$se) && record$se > 0 &&
    is.finite(record$p) && record$p >= 0 && record$p <= 1
  record
}
