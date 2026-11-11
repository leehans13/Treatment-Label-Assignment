# Resume completed work only when the analysis inputs and implementation match.
make_manifest <- function(inputs, units, scripts, extra = list()) {
  info <- file.info(inputs)
  list(paths = normalizePath(inputs), sizes = info$size, modified = as.numeric(info$mtime),
       units = units, extra = extra, scripts = unname(tools::md5sum(scripts)),
       packages = vapply(c("nebula", "fixest", "glmGamPoi", "edgeR", "limma", "Matrix", "SingleCellExperiment"),
                         function(p) as.character(utils::packageVersion(p)), character(1)),
       r_version = as.character(getRversion()), nebula_ncore = getOption("degpaper.nebula.ncore", 2L))
}
checkpoint_map <- function(units, output, manifest, compute, valid, legacy = NULL) {
  dir <- file.path(output, "checkpoints")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, "manifest.rds")
  if (file.exists(path)) {
    saved_manifest <- readRDS(path)
    if (!identical(saved_manifest, manifest) && (is.null(legacy) || !legacy(saved_manifest, manifest)))
      stop("Checkpoint inputs or software have changed; use a new output directory")
  } else {
    if (length(list.files(dir, pattern = "^target_[0-9]+\\.rds$")))
      stop("Checkpoints exist without a manifest; use a new output directory")
    saveRDS(manifest, path)
  }
  lapply(seq_along(units), function(i) {
    path <- file.path(dir, sprintf("target_%03d.rds", i))
    if (file.exists(path)) {
      saved <- readRDS(path)
      if (!identical(saved$target, units[i]) || !valid(saved$result)) stop("Invalid checkpoint: ", path)
      message("Using checkpoint: ", saved$target)
      return(saved$result)
    }
    value <- compute(units[i])
    if (!valid(value)) stop("Invalid result for: ", units[i])
    tmp <- tempfile(pattern = "checkpoint_", tmpdir = dir)
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(list(target = units[i], result = value), tmp)
    if (!file.rename(tmp, path)) stop("Could not save checkpoint: ", path)
    value
  })
}
