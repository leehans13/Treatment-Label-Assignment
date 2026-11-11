# Recreate the active SCT count assay and relevant metadata from h5Seurat.
read_schmidt <- function(path) {
  f <- hdf5r::H5File$new(path, mode = "r")
  on.exit(f$close_all())
  assay <- hdf5r::h5attr(f, "active.assay")
  if (assay != "SCT") stop("Expected the SCT assay; found ", assay)
  group <- f[[paste0("assays/", assay, "/counts")]]
  dims <- as.integer(hdf5r::h5attr(group, "dims"))
  cells <- f[["cell.names"]][]
  genes <- f[[paste0("assays/", assay, "/features")]][]
  mat <- methods::new("dgCMatrix", Dim = dims,
    p = as.integer(group[["indptr"]][]), i = as.integer(group[["indices"]][]),
    x = as.numeric(group[["data"]][]), Dimnames = list(genes, cells))
  g <- f[["meta.data/gene"]]
  gene <- g[["levels"]][][g[["values"]][]]
  meta <- data.frame(gene = gene, donor = f[["meta.data/donor"]][], row.names = cells)
  stopifnot(nrow(meta) == ncol(mat), !anyNA(meta$gene), !anyNA(meta$donor))
  list(counts = mat, metadata = meta)
}

schmidt_genes <- function(counts, metadata, cytokine_path) {
  cytokines <- sort(unique(readr::read_tsv(cytokine_path, col_names = FALSE, show_col_types = FALSE)$X3))
  granzymes <- c("GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1", "GNLY", "THBS1")
  test <- c(cytokines[cytokines %in% rownames(counts)], granzymes[granzymes %in% rownames(counts)])
  target <- unique(as.character(metadata$gene[metadata$gene != "NO-TARGET"]))
  # Select normalization genes using SCT counts across all cells.
  genes <- unique(c(test, rownames(counts)[Matrix::rowMeans(counts) > 1]))
  list(test = test, targets = target, analysis = genes)
}

analyze_schmidt_target <- function(counts, metadata, target, genes) {
  selected <- metadata$gene %in% c("NO-TARGET", target)
  meta <- metadata[selected, c("gene", "donor"), drop = FALSE]
  names(meta) <- c("tx_cell", "donor_id")
  meta$tx_cell <- as.integer(meta$tx_cell != "NO-TARGET")
  meta <- meta[order(meta$donor_id), , drop = FALSE]
  meta$id <- meta$donor_id
  cnt <- counts[genes, rownames(meta), drop = FALSE]
  cnt <- cnt[Matrix::rowMeans(cnt) > 0.01, , drop = FALSE]
  message(target, ": ", nrow(cnt), " genes, ", ncol(cnt), " cells")
  unname(run_models(cnt, meta, "donor_id"))
}

summarize_schmidt <- function(result, targets, test_genes, output) {
  stopifnot(length(result) == length(targets))
  labels <- model_labels
  rows <- lapply(seq_along(result), function(i) {
    methods <- lapply(result[[i]], function(x) { x <- as.data.frame(x); x$gene <- rownames(x); x })
    names(methods) <- labels
    data.table::rbindlist(methods, idcol = "method")
  })
  names(rows) <- targets
  df <- data.table::rbindlist(rows, idcol = "perturb_gene") |>
    dplyr::filter(gene %in% test_genes) |>
    dplyr::mutate(gene_perturb = paste(gene, perturb_gene, sep = "_")) |>
    dplyr::select(gene_perturb, method, `Pr(>|t|)`) |>
    tidyr::pivot_wider(names_from = method, values_from = `Pr(>|t|)`)
  write_summaries(df, 0.05 / 4807, output, "schmidt")
}
