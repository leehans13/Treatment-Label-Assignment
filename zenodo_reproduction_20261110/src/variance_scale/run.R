#!/usr/bin/env Rscript
# Optional arguments: output directory, replicates, bootstrap draws, mode, worker count.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 5L) stop('Usage: Rscript src/variance_scale/run.R [output] [replicates] [bootstrap] [full|smoke] [workers]')
output <- if (length(args)) args[1] else 'data/variance_scale_100'
reps <- if (length(args) >= 2) as.integer(args[2]) else 100L
boot <- if (length(args) >= 3) as.integer(args[3]) else 2000L
mode <- if (length(args) >= 4) args[4] else 'full'
workers <- if (length(args) >= 5) as.integer(args[5]) else 4L
if (is.na(reps) || reps < 2 || is.na(boot) || boot < 2 || !mode %in% c('full','smoke') || is.na(workers) || workers < 1) stop('Invalid run parameters')
source('src/common/checkpoints.R'); source('src/variance_scale/analysis.R')
source('src/variance_scale/summarize.R')
options(degpaper.nebula.ncore = 1L)
grid <- expand.grid(scenario = 1:2, n = if (mode == 'full') c(10L,20L,50L,100L) else c(10L,20L),
                    m = if (mode == 'full') c(20L,50L,100L,200L) else 20L)
set.seed(20261108L)
seeds <- sample.int(.Machine$integer.max, nrow(grid) * reps)
scripts <- c('src/variance_scale/analysis.R', 'src/variance_scale/run.R')
dir.create(output, recursive = TRUE, showWarnings = FALSE)
# Each setting owns its own checkpoint directory, so workers never write the same file.
run_setting <- function(setting) {
  g <- grid[setting, ]
  units <- sprintf('%03d-%03d', setting, seq_len(ceiling(reps / 25)))
  checkpoint_dir <- file.path(output, 'settings', sprintf('%03d', setting))
  manifest <- make_manifest(character(), units, scripts,
    extra = list(grid = grid, reps = reps, seed = 20261108L, rng_kind = RNGkind(), nebula = as.character(packageVersion('nebula'))))
  checkpoint_map(units, checkpoint_dir, manifest, compute = function(unit) {
    block <- as.integer(strsplit(unit, '-', fixed = TRUE)[[1]][2])
    idx <- seq.int((block - 1L) * 25L + 1L, min(block * 25L, reps))
    message('Scenario ', g$scenario, ', n=', g$n, ', m=', g$m, ', replicates ', min(idx), '-', max(idx))
    do.call(rbind, lapply(idx, function(i) {
      result <- simulate_fit(g$scenario, g$n, g$m, seeds[(setting - 1L) * reps + i])
      result$replicate <- i; result
    }))
  }, valid = function(x) is.data.frame(x) && nrow(x) > 0L)
}
# Per-replicate seeds make results independent of worker count and scheduling.
if (.Platform$OS.type == 'windows' && workers > 1L) stop('Parallel workers require a Unix-like system')
results <- parallel::mclapply(seq_len(nrow(grid)), run_setting, mc.cores = min(workers, nrow(grid)), mc.preschedule = FALSE, mc.set.seed = FALSE)
failed <- vapply(results, inherits, logical(1), what = 'try-error')
if (any(failed)) stop('Settings failed: ', paste(which(failed), collapse = ', '), '; see log. Completed settings remain checkpointed.')
data <- do.call(rbind, unlist(results, recursive = FALSE))
saveRDS(data, file.path(output, 'replicates.rds'))
write.csv(data, file.path(output, 'replicates.csv'), row.names = FALSE)
summarize_variance(data, output, boot)
