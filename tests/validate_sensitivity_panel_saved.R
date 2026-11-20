# Independent author-run saved-panel verification. No fitting or q/P adjustment.
validate_sensitivity_panel_saved <- function(root,run,export) {
  source(file.path(root,"tests/release_validation_config.R"),local=environment())
  config<-release_publication_config(root)
  b<-publication_data(run,config);sens<-publication_sensitivity_pairs(b,run)
  specs<-list(A=c("BMI23","species"),B=c("CLR","species"),C=c("CLR","pathways"))
  index<-utils::read.csv(file.path(export,"figure_index.csv"),stringsAsFactors=FALSE)
  dimensions<-yaml::read_yaml(file.path(root,"config/publication_panel_config.yml"))$sensitivity_scatter
  all_files<-list.files(file.path(export,"figures"),recursive=TRUE,full.names=FALSE)
  assert_true(!any(grepl("(^|/)SFig5|diversity_diagnostics",all_files)),"Obsolete supplementary picture in fresh export")
  read<-function(n)utils::read.csv(file.path(export,"plot_data",n),check.names=FALSE,stringsAsFactors=FALSE)
  for(letter in names(specs)) {
    spec<-specs[[letter]];id<-paste0("SFig4",letter,"_sensitivity_",paste(spec,collapse="_"))
    row<-index[index$figure_id==id,,drop=FALSE]
    assert_true(nrow(row)==1L&&row$width_mm==dimensions$width_mm&&row$height_mm==dimensions$height_mm,"Sensitivity export dimensions differ from current settings")
    expected<-sens$pairs[sens$pairs$method==spec[1]&sens$pairs$domain==spec[2],,drop=FALSE]
    saved<-read(paste0(id,".csv"))
    assert_true(nrow(saved)==nrow(expected)&&!anyDuplicated(saved[c("comparison","feature")]),"Matched feature set changed")
    # Validate original fields/order as well as the added display variables.
    for(n in names(expected))assert_true(isTRUE(all.equal(saved[[n]],expected[[n]],tolerance=1e-12,check.attributes=FALSE)),paste("Sensitivity plot changed source field:",id,n))
    assert_true(identical(saved$direction_class,ifelse(expected$direction_agrees,"Concordant","Discordant")),"Concordance colours mismatch")
    shared<-!is.na(expected$q_main)&!is.na(expected$q_sensitivity)&expected$q_main<.05&expected$q_sensitivity<.05
    assert_true(identical(saved$shared_fdr,shared),"Shared FDR definition changed")
    assert_true(identical(saved$shared_status,ifelse(shared,"Shared FDR-significant","Not shared FDR-significant")),"Shared alpha class mismatch")
    annotations<-read(paste0(id,"_facet_summary.csv"))
    original<-sens$summary[sens$summary$method==spec[1]&sens$summary$domain==spec[2],,drop=FALSE]
    assert_true(identical(annotations$comparison,original$comparison),"Facet summary order changed")
    for(n in names(original))assert_true(isTRUE(all.equal(annotations[[n]],original[[n]],tolerance=1e-12,check.attributes=FALSE)),paste("Descriptive facet summary differs:",n))
    assert_true(all(grepl("n = ",annotations$label,fixed=TRUE))&&all(grepl("Concordant = ",annotations$label,fixed=TRUE))&&
      all(grepl("r = ",annotations$label,fixed=TRUE))&&all(grepl("Shared FDR n = ",annotations$label,fixed=TRUE)),"Required facet annotation missing")
    palette<-read(paste0(id,"_direction_palette.csv"));alpha<-read(paste0(id,"_shared_alpha.csv"))
    assert_true(identical(palette$direction_class,c("Concordant","Discordant"))&&identical(palette$colour,c("#2166AC","#F7AA58")),"Sensitivity direction palette differs")
    assert_true(identical(alpha$shared_status,c("Not shared FDR-significant","Shared FDR-significant"))&&identical(alpha$alpha,c(.18,.95)),"Sensitivity point opacity differs")
  }
  invisible(TRUE)
}
