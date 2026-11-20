annotate_pathway_features <- function(features, root) {
  base <- sub("[|].*$", "", features)
  id <- sub(":.*$", "", base)
  description <- trimws(sub("^[^:]+:?", "", base))
  has_genus <- grepl("g__", features, fixed = TRUE)
  has_species <- grepl("s__", features, fixed = TRUE)
  genus <- ifelse(has_genus, sub("[.|;].*$", "", sub("^.*g__", "", features)), NA_character_)
  species <- ifelse(has_species, sub("[|;].*$", "", sub("^.*s__", "", features)), NA_character_)
  contributor <- ifelse(has_species, species, genus)
  dictionary <- utils::read.csv(file.path(root, "config/pathway_function_dictionary.csv"), check.names = FALSE)
  category <- rep("Other/Unclassified", length(features))
  for (i in seq_len(nrow(dictionary))) {
    hit <- category == "Other/Unclassified" &
      grepl(dictionary$pattern[i], gsub("[._]+", " ", description), ignore.case = TRUE)
    category[hit] <- dictionary$category[i]
  }
  data.frame(pathway_feature = features, base_pathway_ID = id, pathway_description = description,
    taxonomic_level = ifelse(has_species, "species", ifelse(has_genus, "genus", "unstratified")),
    genus = genus, taxonomic_contributor = gsub("_", " ", contributor),
    display_name = base, functional_category = category,
    annotation_source = "Provider HUMAnN taxonomic string; fixed keyword classification dictionary",
    classification_status = ifelse(category == "Other/Unclassified", "Unclassified", "Dictionary match"))
}

select_pathway_candidates <- function(x, config) {
  stopifnot(all(c("coef", "qval", "pathway") %in% names(x)))
  annotation <- annotate_pathway_features(x$pathway, config$project_root)
  d <- cbind(x, annotation[, setdiff(names(annotation), "pathway_feature"), drop = FALSE])
  d$visualization_candidate <- abs(d$coef) > 1 & d$qval < .01
  d$selection_rank <- NA_integer_
  d$selected_for_panel <- FALSE
  d$genus_panel_rank <- NA_integer_
  d$species_panel_rank <- NA_integer_
  d$selected_for_genus_panel <- FALSE
  d$selected_for_species_panel <- FALSE
  for (rank in c("genus", "species")) {
    # A species-stratified HUMAnN feature also has a genus assignment. Figure
    # 5B summarizes that assignment; it is not a separately fitted genus model.
    eligible <- if (rank == "genus") !is.na(d$genus) & nzchar(d$genus) else d$taxonomic_level == "species"
    ix <- which(d$visualization_candidate & eligible)
    ix <- ix[order(d$qval[ix], -abs(d$coef[ix]), d$feature[ix])]
    d$selection_rank[ix] <- seq_along(ix)
    take <- head(ix, config$publication$pathway_candidate_limit[[rank]])
    d$selected_for_panel[take] <- TRUE
    d[[paste0(rank, "_panel_rank")]][ix] <- seq_along(ix)
    d[[paste0("selected_for_", rank, "_panel")]][take] <- TRUE
  }
  d
}
