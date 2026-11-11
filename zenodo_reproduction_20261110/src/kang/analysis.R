# Analysis logic adapted from perturb/kang.ipynb.
prepare_kang <- function(obj) {
  stopifnot(inherits(obj, "SingleCellExperiment"))
  cd <- as.data.frame(SummarizedExperiment::colData(obj))
  # Map ExperimentHub metadata to the analysis column names.
  if (!all(c("group_id", "patient_id", "cluster_id") %in% names(cd))) {
    required <- c("stim", "ind", "cell")
    if (!all(required %in% names(cd))) stop("Input lacks the required Kang cell metadata")
    cd$group_id <- cd$stim
    cd$patient_id <- cd$ind
    cd$cluster_id <- cd$cell
  }
  if (anyNA(cd$group_id) || anyNA(cd$patient_id)) stop("Missing treatment or donor metadata")
  if (!all(as.character(cd$group_id) %in% c("ctrl", "stim"))) stop("Unexpected treatment labels")
  if (anyNA(cd$cluster_id)) {
    message("Dropping ", sum(is.na(cd$cluster_id)), " cells without a cell-type label")
    keep <- !is.na(cd$cluster_id)
    obj <- obj[, keep]
    cd <- cd[keep, , drop = FALSE]
  }
  if (!"counts" %in% SummarizedExperiment::assayNames(obj)) stop("Input has no counts assay")
  SummarizedExperiment::colData(obj)$group_id <- cd$group_id
  SummarizedExperiment::colData(obj)$patient_id <- cd$patient_id
  SummarizedExperiment::colData(obj)$cluster_id <- cd$cluster_id
  obj
}

analyze_kang <- function(obj, map_units = function(units, compute) lapply(units, compute)) {
  cd <- as.data.frame(SummarizedExperiment::colData(obj))
  cell_types <- unique(as.character(cd$cluster_id))
  results <- map_units(cell_types, function(cell_type) {
    message("Analyzing ", cell_type)
    selected <- which(as.character(cd$cluster_id) == cell_type)
    meta <- cd[selected, c("group_id", "patient_id", "cluster_id"), drop = FALSE]
    names(meta) <- c("tx_cell", "donor_id", "cell_type")
    meta$tx_cell <- as.integer(as.character(meta$tx_cell) != "ctrl")
    meta <- meta[order(meta$donor_id), , drop = FALSE]
    meta$id <- meta$donor_id
    cnt <- SummarizedExperiment::assay(obj, "counts")[, rownames(meta), drop = FALSE]
    cnt <- cnt[Matrix::rowMeans(cnt) > 0.01, , drop = FALSE]
    message("  ", nrow(cnt), " genes, ", ncol(cnt), " cells")
    out <- lapply(run_models(cnt, meta, "id"), function(value) {
      value <- as.data.frame(value)
      value$gene <- rownames(value)
      value
    })
    data.table::rbindlist(out, idcol = "method")
  })
  names(results) <- cell_types
  data.table::rbindlist(results, idcol = "celltype")
}

summarize_kang <- function(result, output_dir) {
  # Use the same method label in the output summaries.
  result$method[result$method == "robust GLM (cell)"] <- "fixest (cell)"
  df <- result |>
    dplyr::mutate(gene_celltype = paste(gene, celltype, sep = "_")) |>
    dplyr::select(gene_celltype, method, `Pr(>|t|)`) |>
    tidyr::pivot_wider(names_from = method, values_from = `Pr(>|t|)`)
  write_summaries(df, 0.05 / 6257, output_dir, "kang")
}
