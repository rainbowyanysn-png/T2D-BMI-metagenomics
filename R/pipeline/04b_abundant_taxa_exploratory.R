run_abundant_taxa_exploratory <- function(data, run, config) {
  prepared <- abundant_taxa_prepare(data, run, config)
  results <- abundant_taxa_fit(prepared, config)
  for (name in names(results)) safe_write_csv(results[[name]],
    path_in_run(run, "tables/04b_abundant_taxa", paste0(name, ".csv")))
  invisible(results)
}
