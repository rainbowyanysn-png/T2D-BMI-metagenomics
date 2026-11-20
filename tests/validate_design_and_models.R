validate_design <- function(metadata, fixed, references) {
  assert_columns(metadata, c("sample_id", fixed), "Model metadata")
  assert_true(!anyNA(metadata[, fixed, drop = FALSE]), "Missing covariates before model fitting.")
  x <- as.data.frame(metadata)
  for (reference in references) {
    parts <- strsplit(reference, ",", fixed = TRUE)[[1]]
    x[[parts[1]]] <- stats::relevel(factor(x[[parts[1]]]), ref = parts[2])
  }
  design <- stats::model.matrix(stats::reformulate(fixed), data = x)
  assert_true(qr(design)$rank == ncol(design), "Rank-deficient design matrix.")
  list(matrix = design, metadata = x, references = references)
}
