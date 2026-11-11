# Summaries use whole-replicate bootstrap; binomial intervals describe Monte Carlo uncertainty.
summarize_variance <- function(data, output, B = 2000L) {
  set.seed(20261109L)
  groups <- split(data, interaction(data$scenario, data$n, data$m, drop = TRUE))
  interval <- function(k, n) if (n > 0) unname(stats::binom.test(k, n)$conf.int) else c(NA_real_, NA_real_)
  quants <- function(x) if (sum(is.finite(x)) >= 2) unname(quantile(x[is.finite(x)], c(.025,.975))) else c(NA_real_, NA_real_)
  summaries <- lapply(groups, function(d) {
    ok <- d[d$valid, , drop = FALSE]; count <- nrow(ok); total <- nrow(d)
    fail_ci <- interval(total-count, total)
    boundary_known <- is.finite(d$subject_variance)
    boundary_n <- sum(d$boundary[boundary_known], na.rm = TRUE)
    boundary_ci <- interval(boundary_n, sum(boundary_known))
    rejected <- sum(ok$p < .05); covered <- sum(abs(ok$estimate) <= qnorm(.975)*ok$se)
    reject_ci <- interval(rejected,count); cover_ci <- interval(covered,count)
    empirical <- if (count >= 2) var(ok$estimate) else NA_real_
    bs <- if (count >= 2) replicate(B, {
      i <- sample.int(count, count, replace = TRUE)
      var(ok$estimate[i])
    }) else rep(NA_real_, B)
    ev_ci <- quants(bs)
    data.frame(scenario=d$scenario[1], n=d$n[1], m=d$m[1], total=total, valid=count,
      failures=total-count, failure_rate=(total-count)/total, failure_low=fail_ci[1], failure_high=fail_ci[2],
      boundary_count=boundary_n, boundary_denominator=sum(boundary_known), boundary_low=boundary_ci[1], boundary_high=boundary_ci[2],
      empirical=empirical, empirical_low=ev_ci[1], empirical_high=ev_ci[2],
      rejection=if(count) rejected/count else NA_real_, rejection_low=reject_ci[1], rejection_high=reject_ci[2],
      coverage=if(count) covered/count else NA_real_, coverage_low=cover_ci[1], coverage_high=cover_ci[2])
  })
  summary <- do.call(rbind, summaries); rownames(summary) <- NULL
  summary <- summary[order(summary$scenario,summary$n,summary$m), ]
  write.csv(summary, file.path(output,'summary.csv'), row.names=FALSE)
  reductions <- function(by, varying) {
    groups <- split(summary, interaction(summary$scenario, summary[[by]], drop=TRUE))
    rows <- lapply(groups, function(z) {
      z <- z[order(z[[varying]]), ]; if (nrow(z)<2) return(NULL)
      data.frame(scenario=z$scenario[-1], fixed_variable=by, fixed_value=z[[by]][-1], varying_variable=varying,
        from=z[[varying]][-nrow(z)], to=z[[varying]][-1],
        relative_reduction=1-z$empirical[-1]/z$empirical[-nrow(z)])
    })
    do.call(rbind,rows)
  }
  red <- rbind(reductions('n','m'),reductions('m','n'))
  if (!is.null(red)) write.csv(red,file.path(output,'variance_reductions.csv'),row.names=FALSE)
  invisible(summary)
}
