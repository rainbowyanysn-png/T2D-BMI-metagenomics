load_analysis_support <- function(root,scope) {
 for(n in c("assertions","model_helpers","clinical_correlations","pooled_pathway_maaslin","group_display","source_resources","data_import","diversity_models","abundant_taxa_models","publication_display_data","medication_audit","metadata_audit","checkpoints","species_loco_helpers","figure_exports","pipeline_runner"))
  source(file.path(root,"R/functions",paste0(n,".R")),local=scope)
}
