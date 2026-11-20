final_beta_sources <- function(run) {
  read<-function(n)pub_read(run,paste0("tables/04c_raw_Aitchison/BMI25_raw_Aitchison_",n,".csv"))
  disease<-read("PERMDISP_disease");groups<-read("PERMDISP_four_groups")
  extract<-function(x,group,expected) {
    z<-x[x$term=="Groups",,drop=FALSE]
    assert_true(nrow(z)==1L&&abs(z$F-expected)<1e-12&&z[["Pr(>F)"]]==.001&&z$type=="Spatial median"&&z$n_samples==1182L&&z$n_permutations==999L,paste("Formal raw-Aitchison PERMDISP differs:",group))
    z$comparison<-group;z$inference_role<-"Formal; uncorrected CZM/CLR Aitchison";z
  }
  dispersion<-pub_bind(list(extract(disease,"T2D vs Control",22.4781054358386),extract(groups,"BMI25 four groups",6.65286199089215)))
  marginal<-pub_bind(lapply(c("total","four_groups"),function(s){z<-read(paste0("PERMANOVA_",s,"_blocked"));z$comparison<-if(s=="total")"T2D vs Control"else"BMI25 four groups";z$inference_role<-"Formal marginal PERMANOVA; uncorrected Aitchison";z}))
  list(dispersion=dispersion,marginal=marginal,cohort=read("PERMANOVA_by_cohort"))
}

final_beta_caption <- function(run) {
  x<-final_beta_sources(run)$dispersion
  sprintf("Formal beta-diversity inference uses uncorrected CZM/CLR Aitchison distances. Spatial-median PERMDISP: T2D vs Control F = %.3f, P = %.3f; BMI25 four groups F = %.3f, P = %.3f. Dispersion differed, so apparent ordination separation cannot be attributed exclusively to centroid differences. Full model and permutation details are provided in Supplementary Table S3.",x$F[1],x[["Pr(>F)"]][1],x$F[2],x[["Pr(>F)"]][2])
}

final_beta_sheets <- function(run) {
  beta<-final_beta_sources(run)
  module<-function(m,d,p)pub_read_module(run,m,paste0("tables/",d),p)
  combat<-module("primary_diversity","04_primary_diversity","batch_effect_permanova|beta_permanova_global")
  combat$inference_role<-"Descriptive only; ComBat-corrected, not formal beta-diversity inference"
  list(Cohort_PERMANOVA=beta$cohort,Marginal_PERMANOVA=beta$marginal,PERMDISP=beta$dispersion,
    Interactions=module("diversity_interactions","05_diversity_interactions","beta_(BMI|interaction)"),
    BMI23=module("bmi23_diversity","07_BMI23_diversity","Beta_Aitchison_PERM"),Descriptive_ComBat=combat,
    Raw_Aitchison=data.frame(n_samples=nrow(pub_read(run,"tables/04c_raw_Aitchison/BMI25_raw_Aitchison_sample_manifest.csv")),
      n_features=nrow(pub_read(run,"tables/04c_raw_Aitchison/BMI25_raw_Aitchison_feature_manifest.csv")),distance="Uncorrected CZM/CLR Aitchison",zero_handling="CZM; prespecified 80% zero deletion rule"),
    Notes=data.frame(Item=c("Formal inference","PERMDISP","ComBat","Sample/feature handling","Interpretation"),
      Definition=c("Uncorrected CZM/CLR Aitchison; marginal PERMANOVA and spatial-median PERMDISP",final_beta_caption(run),
        "Descriptive visualization only; ComBat results are separated from formal tests","1182 samples and 138 features after the prespecified CZM rule; no threshold changed",
        "Significant dispersion can contribute to apparent separation; no confirmed effect modification is claimed")))
}

final_loco_sheets <- function(run) {
  read<-function(n)pub_read(run,paste0("tables/02a_species_LOCO/",n,".csv"))
  list(LOCO_primary_significant=read("LOCO_primary_significant"),LOCO_feature_summary=read("LOCO_feature_summary"),
    LOCO_all_results=read("LOCO_all_results"),LOCO_model_summary=read("LOCO_model_summary"),
    LOCO_exclusions=read("LOCO_prevalence_exclusions"),LOCO_data_dictionary=read("LOCO_data_dictionary"))
}

final_manuscript_beta_text <- function(run,export) {
  text<-final_beta_caption(run)
  directory<-file.path(export,"manuscript_updates");dir.create(directory,recursive=TRUE,showWarnings=FALSE)
  safe_write_lines(c("Results paragraph for author review",text),file.path(directory,"Beta_diversity_Results.txt"))
  safe_write_lines(c("Response Letter paragraph for author review",
    paste("We distinguish formal uncorrected Aitchison inference from descriptive ComBat ordination.",text,
      "Figure 2 captions and Supplementary Table S3 report the same current-run dispersion statistics.")),file.path(directory,"Response_Letter_beta_diversity.txt"))
  safe_write_lines(c("Author action: insert the reviewed paragraphs into the manuscript/response letter. No external document was modified.",
    "Figure2F and supplementary beta diagnostic plots are excluded; full results remain in S3."),file.path(directory,"README.txt"))
}
