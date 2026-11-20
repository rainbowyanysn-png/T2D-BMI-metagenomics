validate_inputs <- function(data, expected) {
  assert_equal(nrow(data$metadata), expected$inputs$n_samples, "validation: n samples")
  assert_equal(nrow(data$species), expected$inputs$n_species, "validation: n species")
  assert_equal(nrow(data$pathways), expected$inputs$n_pathways, "validation: n pathways")
  invisible(TRUE)
}
