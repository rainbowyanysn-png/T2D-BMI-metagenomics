publication_table_plan <- function(b, run, config) {
  read <- function(p) pub_read(run, paste0("tables/", p))
  module <- function(m, d, pattern = "[.]csv$") pub_read_module(run, m, paste0("tables/", d), pattern)
  note <- "BMI25 groups: NC, non-overweight control; OC, overweight control; NT2D, non-overweight T2D; OT2D, overweight T2D. BMI <25 versus >=25 kg/m²."
  definitions <- function(model, q) data.frame(model = model, formal_q_definition = q)
  plan <- list()
  put <- function(id, title, sheets, footnote, main = FALSE) {
    plan[[id]] <<- list(title = title, sheets = sheets, note = footnote, main = main)
  }
  md <- b$md; display <- unname(group_display[md$BMI25_group])
  table1 <- lapply(c("Overall", "NC", "OC", "NT2D", "OT2D"), function(g) {
    d <- md[if (g == "Overall") rep(TRUE, nrow(md)) else display == g, ]
    values <- c(nrow(d), sprintf("%.2f (%.2f)", mean(d$age), stats::sd(d$age)),
      sprintf("%d (%.1f%%)", sum(d$sex == "female"), 100*mean(d$sex == "female")),
      sprintf("%.2f (%.2f)", mean(d$BMI), stats::sd(d$BMI)),
      as.character(table(factor(d$study_name, levels = unlist(config$cohorts$names)))))
    values
  })
  table1 <- data.frame(Characteristic = c("N", "Age, mean (SD)", "Female, n (%)", "BMI, mean (SD)", unlist(config$cohorts$names)),
    setNames(as.data.frame(table1), c("Overall", "NC", "OC", "NT2D", "OT2D")), check.names = FALSE)
  cohort_heading <- table1[1, , drop = FALSE]; cohort_heading[1, ] <- ""; cohort_heading$Characteristic <- "Study cohort, n"
  table1 <- rbind(table1[1:4, ], cohort_heading, table1[-(1:4), ])
  put("Main_Table_1", "Participant characteristics", list(Characteristics = table1),
    paste(note, "Missing age values were imputed using the within-cohort median before analysis. Values are descriptive and are presented as mean (SD) or n (%), as appropriate. MetaCardis controls are non-T2D controls."), TRUE)

  summary <- list()
  record <- function(x, level, comparison, beta, q, model, transform, family) {
    sig <- !is.na(x[[q]]) & x[[q]] < .05
    parts <- strsplit(comparison, " vs ", fixed = TRUE)[[1]]
    data.frame(analysis_level = level, comparison = comparison,
      reference_group = if (length(parts)==2) parts[2] else "See model definition",
      comparison_group = parts[1], n_tested = nrow(x), n_significant = sum(sig),
      positive_significant = sum(sig & x[[beta]] > 0), negative_significant = sum(sig & x[[beta]] < 0),
      model = model, transformation = transform, formal_q_column = q, q_family = family,
      q_family_size = if ("q_family_size" %in% names(x)) unique(x$q_family_size)[1] else nrow(x))
  }
  for (i in 1:5) {
    label <- b$spec$comparison[i]
    s <- b$species[b$species$comparison==label, ]; p <- b$reml[b$reml$comparison==label, ]
    summary[[paste0("s",i)]] <- record(s,"BMI25 species",label,"coef","qval","MaAsLin2, contrast-specific frozen covariates","TSS + LOG; standardize TRUE","Full-model feature-by-covariate coefficients")
    summary[[paste0("p",i)]] <- record(p,"BMI25 pathway features",label,"combined_log2FC","padj","Cohort DESeq2 + REML","Count-like input, poscounts, local dispersion","All successful pooled features within contrast")
  }
  summary$pooled <- record(b$pooled,"Pooled pathway features","T2D vs Control","coef","qval",
    "MaAsLin2 disease + BMI + age + sex + study_name; all fixed effects","TSS + LOG; standardize TRUE","Full-model feature-by-covariate coefficients")
  summary$clinical <- record(b$clinical,"Exploratory clinical correlations","BMI, age and sex","spearman_rho","q_bh_60",
    "Spearman, exact FALSE; pairwise complete; no covariate adjustment","Relative abundance without additional transformation","All 60 displayed tests")
  for (n in unique(paste(b$interactions$BMI_definition,b$interactions$domain))) {
    x <- b$interactions[paste(b$interactions$BMI_definition,b$interactions$domain)==n, ]
    summary[[n]] <- record(x,paste("Formal interaction",x$domain[1]),x$BMI_definition[1],"coef","formal_q",
      "Frozen T2D x BMI model","TSS + LOG; standardize FALSE","Target interaction term within BMI definition and feature domain")
  }
  for (term in c("T2D","BMI_per5")) summary[[term]] <- record(b$additive,"Additive species",term,paste0("coef_",term),
    paste0("q_within_exposure_BH_",term),"disease + BMI_per5 + age + sex + study_name","TSS + LOG; standardize FALSE","Exposure-specific BH across 200 species")
  put("Analysis_Summary_Table","Analysis models and feature families",list(Analysis_summary=pub_bind(summary)),
    paste(note,"Pooled individual pathway and cohort REML analyses address complementary association summaries and have separate testing families."),TRUE)

  # Each closure is evaluated separately: a missing source CSV does not prevent
  # export of independent workbooks. The overall exporter still returns failure.
  plan$Supplementary_Table_S1 <- function() {
    source <- pub_source_review(run$root,config,"cohort",run=run)
    cells <- pub_bind(lapply(c("BMI25","BMI23"),function(threshold) {
      z <- as.data.frame(table(cohort=md$study_name,group=md[[paste0(threshold,"_group")]]))
      names(z)[3] <- "n_samples"; z$BMI_definition <- threshold; z
    }))
    source$BMI25_rule <- "BMI <25 versus >=25 kg/m²"
    source$BMI23_rule <- "QinJ_2012/YuJ_2015: BMI <23 versus >=23; other cohorts retain 25"
    list(title="Cohort sources and sample flow",sheets=list(Cohort_sources=source,
      Sample_flow=read("01b_sample_flow/BMI25_BMI23_pre_QC_flow.csv"),Four_group_counts=cells,
      Final_sample_IDs=md[,c("sample_id","study_name","disease","BMI25_group","BMI23_group")]),
      note=paste(note,"Source evidence grades and individual-level mapping limitations are preserved as documented. No sample is excluded in response to annotation limitations."),main=FALSE)
  }
  put("Supplementary_Table_S2","Alpha diversity models and comparisons",list(
    BMI25_models=module("alpha_details","04d_alpha_details","fixed_coefficients"),
    Variance_components=module("alpha_details","04d_alpha_details","variance_components"),
    Summary_and_ICC=read("04_primary_diversity/BMI25_alpha_LMM_summary.csv"),
    Pairwise=publication_alpha_pairwise(run),
    Interactions=module("diversity_interactions","05_diversity_interactions","^alpha_"),
    BMI23=publication_bmi23_alpha(run), Notes=publication_alpha_notes()),
    paste(note,"BMI25 pairwise estimates are adjusted means, group1 minus group2. p.value is Tukey-adjusted; p.adj applies BH to the six Tukey-adjusted values within each metric and determines pairwise significance. Categorical BMI25 and continuous BMI interactions use separate BH families of two metrics; see Notes for BMI23 models and adjustments."))
  put("Supplementary_Table_S3","Supplementary Table S3. Beta-diversity analyses and dispersion diagnostics",final_beta_sheets(run),final_beta_caption(run))
  species_columns <- intersect(c("feature","comparison","coef","stderr","pval","qval","N","N.not.zero",
    "comparison_direction","reference_group","comparison_group","model_formula","q_family","q_family_size"), names(b$species))
  put("Supplementary_Table_S4","BMI25 species associations",list(All_five_comparisons=b$species[,species_columns,drop=FALSE],
    Model_and_q_definition=definitions("Contrast-specific MaAsLin2; TSS, LOG, standardize TRUE",
      "MaAsLin2 full-model qval; target-term BH is not the formal species criterion")),note)
  plan$Supplementary_Table_S4$sheets <- c(plan$Supplementary_Table_S4$sheets, final_loco_sheets(run))
  plan$Supplementary_Table_S4$note <- paste(note,"LOCO evaluates omission sensitivity, not new feature discovery. Original full-model q values remain formal; target BH is an audit column. Lack of q<0.05 after omission alone does not establish instability. Excluded features have missing estimates and are not counted as direction reversals.")
  clinical <- b$clinical
  clinical$n_missing <- 1186L-clinical$n_complete
  clinical$FDR_significant <- clinical$q_bh_60 < .05
  clinical <- clinical[,intersect(c("feature","display_name","clinical_trait","spearman_rho","p_value","q_bh_60",
    "n_complete","n_missing","FDR_significant"),names(clinical)),drop=FALSE]
  put("Supplementary_Table_S5","Exploratory clinical correlations",list(All_60_tests=clinical,
    Selected_top20_features=b$clinical_selected,Model_and_q_definition=definitions(
      "Spearman, exact FALSE; Female=0 and Male=1; pairwise complete observations",
      "All 60 raw P values corrected together; formal column q_bh_60")),
    "Top 20 features ranked by overall species full-model qval, pval and feature ID. Saved relative abundance without additional transformation; no covariate adjustment. These marginal correlations do not establish independence or causality.")
  pooled <- b$candidates
  pooled_columns <- intersect(c("feature","pathway","base_pathway_ID","pathway_description","taxonomic_contributor",
    "taxonomic_level","coef","stderr","pval","qval","N","N.not.zero","functional_category","annotation_source"),names(pooled))
  pooled_reader <- pooled[,pooled_columns,drop=FALSE]
  put("Supplementary_Table_S6","Pooled multivariable pathway associations",list(All_5945_features=pooled_reader,
    Significant_features=pooled_reader[pooled$qval<.05,],Figure5_candidates=pooled_reader[pooled$visualization_candidate,],
    Model_and_q_definition=definitions("disease + BMI + age + sex + study_name; all fixed effects; TSS + LOG; standardize TRUE",
      "Significance is based on the MaAsLin2 full-model q value across all model coefficients.")),
    "1186 participants; 5945 pathway features. Taxon-stratified contributions reflect inferred functional potential. The display criterion abs(coef)>1 and qval<0.01 does not replace the formal qval<0.05 threshold.")
  display_view <- submission_pathway_view(run$root,pooled)
  put("Supplementary_Table_S7","Taxon-stratified MaAsLin2 pathway features",list(
    Candidates_and_selection=display_view[display_view$visualization_candidate & display_view$taxonomic_level!="unstratified",
      !names(display_view)%in%c("qval_disease_term_BH","audit_target_family_size","selection_rank","genus_panel_rank","species_panel_rank"),drop=FALSE],
    Model_and_q_definition=definitions("Same fitted pooled MaAsLin2 model as Table S6",
      "Full-model qval; all64 eligible species-stratified contributions; six species, four genera and42 pathway IDs; no top-N truncation")),
    "Display flags identify source contributions for Figures5B/5C and SFig3. All eligible species-stratified contributions with full-model qval<0.01 and abs(coef)>1 are displayed; named genus/species/pathway required; UNINTEGRATED/UNMAPPED excluded. Non-displayed candidates remain with explicit reasons. Genus bars aggregate contributions and are not64 separate bars. Assignments indicate inferred functional potential, not measured flux or a causal mechanism.")
  put("Supplementary_Table_S8","Cohort-specific DESeq2 and REML pathway associations",list(
    Pooled_REML_all=b$reml,Pooled_REML_significant=b$reml[b$reml$padj<.05,],Cohort_effects=b$cohort_effects,
    Heterogeneity=b$reml[,c("pathway","comparison","n_studies","Heterogeneity_I2","stderr","ci_lower","ci_upper")],
    Model_and_q_definition=definitions("Cohort DESeq2: counts TRUE, prevalence 20%, poscounts, local dispersion; >=2 valid cohorts; REML",
      "BH-adjusted P value within each complete successful pooled contrast; padj")),note)
  put("Supplementary_Table_S9","Pathway annotation and functional classification",list(
    Feature_annotation=b$reml_annotation,Fixed_dictionary=utils::read.csv(file.path(run$root,"config/pathway_function_dictionary.csv"))),
    "Fixed keyword classification is independent of effect direction. Unmatched features are Other/Unclassified. Features are not assumed to be independent biochemical pathways.")
  interaction_reader <- b$interactions[,!names(b$interactions)%in%c("q_maaslin_full_model","qval","qval_full_model"),drop=FALSE]
  put("Supplementary_Table_S10","Formal T2D by BMI interactions",list(All_six_families=interaction_reader,
    Model_and_q_definition=definitions("Frozen BMI25, continuous BMI and cohort-adapted BMI23 interaction models; standardize FALSE",
      "Target interaction-term BH q value; separate BMI definitions and species/pathway families")),
    "Continuous BMI is centered using all 1186 participants; coefficients are per 1 kg/m² with corresponding per-5 rescaling. Significance uses target interaction-term BH q values within each feature domain and BMI definition.")
  additive_reader<-b$additive[,!grepl("^qval_|^q_joint",names(b$additive)),drop=FALSE]
  additive_long_reader<-b$additive_long[,!names(b$additive_long)%in%c("qval","q_maaslin_full_model","q_joint_disease_BMI_BH"),drop=FALSE]
  put("Supplementary_Table_S11","T2D and BMI additive associations",list(All_200_species=additive_reader,
    Target_coefficients=additive_long_reader,Model_and_q_definition=definitions(
      "disease + BMI_per5 + age + sex + study_name; standardize FALSE; BMI centered on all participants",
      "Exposure-specific BH across 200 species for T2D and BMI separately")),
    "BMI effects are per 5 kg/m². Association classes use separate exposure-specific BH families of 200 species.")
  bmi23s <- module("bmi23_species","07_BMI23_sensitivity","_all[.]csv$")
  bmi23p <- module("bmi23_pathways","07_BMI23_pathways","_all[.]csv$")
  clr <- module("clr","08_CLR_sensitivity","_all[.]csv$")
  sensitivity <- publication_sensitivity_pairs(b,run)
  put("Supplementary_Table_S12","Cohort-adapted BMI23 sensitivity analysis",list(Species=bmi23s,Pathways=bmi23p,
    Formal_interactions=interaction_reader[interaction_reader$BMI_definition=="BMI23",],
    Sample_counts=as.data.frame(table(cohort=md$study_name,group=md$BMI23_group)),
    Concordance=sensitivity$summary[sensitivity$summary$method=="BMI23",]),
    "QinJ_2012 and YuJ_2015 use BMI23; other cohorts retain BMI25. Four-group DESeq2 and diversity factor order follows the source implementation: NC, NT2D, OC, OT2D. Species full-model qval and per-contrast pooled pathway padj remain distinct.")
  put("Supplementary_Table_S13","CZM plus CLR sensitivity analyses",list(Complete_results=clr,
    Concordance=sensitivity$summary[sensitivity$summary$method=="CLR",],
    Paired_species=sensitivity$pairs[sensitivity$pairs$method=="CLR" & sensitivity$pairs$domain=="species",]),
    "CZM replacement, CLR transformation and limma empirical-Bayes models; BH within each target contrast or interaction family. The T2D-internal species model adjusts the same categorical treatment covariate as its primary model; pathway and interaction models do not. BMI23 interaction uses QinJ_2012/YuJ_2015 at 23 and all other cohorts at 25. Coefficient scales and statistical frameworks differ; concordance is descriptive.")
  medication <- list(); medication_summary <- list()
  for (threshold in c("BMI25","BMI23")) {
    u <- read(paste0("09_no_medication_covariate/",threshold,"_all.csv"))
    a <- if (threshold=="BMI25") read("02_primary_BMI25_species/OT2D_vs_NT2D_BMI25_all.csv") else
      read("07_BMI23_sensitivity/OT2D_vs_NT2D_BMI23_all.csv")
    u$BMI_definition <- a$BMI_definition <- threshold
    u$medication_covariate <- FALSE; a$medication_covariate <- TRUE
    u$model <- "Sensitivity model without medication covariate"; a$model <- "Medication-adjusted primary model"
    medication[[paste0(threshold,"_unadjusted")]] <- u
    medication[[paste0(threshold,"_adjusted")]] <- a
    z <- merge(a[,c("feature","coef","qval")],u[,c("feature","coef","qval")],by="feature",suffixes=c("_adjusted","_unadjusted"))
    z$BMI_definition <- threshold; z$direction_agrees <- sign(z$coef_adjusted)==sign(z$coef_unadjusted)
    medication[[paste0(threshold,"_paired")]] <- z
    good <- is.finite(z$coef_adjusted) & is.finite(z$coef_unadjusted)
    medication_summary[[threshold]] <- data.frame(BMI_definition=threshold,n_matched=nrow(z),
      direction_agreement=mean(z$direction_agrees[good]),
      descriptive_coefficient_r=stats::cor(z$coef_adjusted[good],z$coef_unadjusted[good]),
      adjusted_significant=sum(z$qval_adjusted<.05),unadjusted_significant=sum(z$qval_unadjusted<.05),
      adjusted_only=sum(z$qval_adjusted<.05 & z$qval_unadjusted>=.05),
      shared_significant=sum(z$qval_adjusted<.05 & z$qval_unadjusted<.05,na.rm=TRUE))
  }
  internal_medication <- list(title="Medication-adjustment sensitivity",sheets=list(
    Models=pub_bind(medication[!grepl("paired$",names(medication))]),
    Paired_coefficients=pub_bind(medication[grepl("paired$",names(medication))]),
    Concordance=pub_bind(medication_summary)),
    note="Full-model MaAsLin2 q value determines significance. Medication-adjusted models are primary; models without the medication covariate are sensitivity analyses. All paired BMI25 and BMI23 species models use standardize TRUE. This is not a metformin-naive analysis.")
  attr(plan,"internal_medication") <- internal_medication
  plan$Supplementary_Table_S14 <- function() list(title="Clinical and confounder availability",
    sheets=list(Source_and_database=submission_source_view(pub_source_review(run$root,config,"confounder",run=run)),
      Database_completeness=local({z<-read("11_metadata_medication_audit/availability/03_FieldAvailability_ByCohort.csv");
        z[!z$variable_id%in%"historical_antidiabetes_group",,drop=FALSE]})),
    note="Source-level collection and database retention are distinct. Evidence grades and source limitations remain as documented. Categorical antidiabetic treatment was included in the T2D-internal species models. A reliably defined metformin-free analysis was not available. Metformin was not a separate binary covariate. Unknown T2D treatment falls into Other_Anti-diabetic_Drugs under the prespecified classification, not verified nonuse. Pathway meta-analysis and formal interaction models were not medication-adjusted.",main=FALSE)
  plan
}

publication_sensitivity_pairs <- function(b,run) {
  out <- list(); summaries <- list()
  ids <- c("NT2D_vs_NC_BMI23","OT2D_vs_OC_BMI23","OC_vs_NC_BMI23","OT2D_vs_NT2D_BMI23")
  clrids <- c("Total_T2D_vs_control","Lean_T2D_vs_Lean_Control","Overweight_T2D_vs_Overweight_Control",
    "Overweight_vs_Lean_Controls","Overweight_vs_Lean_T2D")
  for (domain in c("species","pathways")) for (i in 1:5) for (method in c("BMI23","CLR")) {
    if (method=="BMI23" && i==1L) next
    main <- if (domain=="species") b$species[b$species$comparison==b$spec$comparison[i],] else b$reml[b$reml$comparison==b$spec$comparison[i],]
    a <- data.frame(feature=make.names(if(domain=="species")main$feature else main$pathway),
      coef_main=if(domain=="species")main$coef else main$combined_log2FC,q_main=if(domain=="species")main$qval else main$padj)
    file <- if (method=="BMI23") {
      if(domain=="species") paste0("tables/07_BMI23_sensitivity/",ids[i-1L],"_all.csv") else
        paste0("tables/07_BMI23_pathways/BMI23_",c("Lean_T2D_vs_Lean_Control","Overweight_T2D_vs_Overweight_Control","Overweight_vs_Lean_Controls","Overweight_vs_Lean_T2D")[i-1L],"_all.csv")
    } else paste0("tables/08_CLR_sensitivity/BMI25_",if(domain=="species")"Taxonomy" else "Pathway","_",clrids[i],"_all.csv")
    z <- pub_read(run,file)
    d <- data.frame(feature=make.names(if(method=="CLR" || domain=="species")z$feature else z$pathway),
      coef_sensitivity=if(method=="CLR")z$coef_CLR else if(domain=="species")z$coef else z$combined_log2FC,
      q_sensitivity=if(method=="CLR")z$qval_CLR_BH else if(domain=="species")z$qval else z$padj)
    stopifnot(!anyDuplicated(a$feature),!anyDuplicated(d$feature))
    paired <- merge(a,d,by="feature"); paired$method<-method;paired$domain<-domain;paired$comparison<-b$spec$comparison[i]
    paired$direction_agrees<-sign(paired$coef_main)==sign(paired$coef_sensitivity)
    key<-paste(method,domain,i);out[[key]]<-paired
    good<-is.finite(paired$coef_main)&is.finite(paired$coef_sensitivity)
    summaries[[key]]<-data.frame(method=method,domain=domain,comparison=b$spec$comparison[i],n_matched=nrow(paired),
      shared_significant=sum(paired$q_main<.05 & paired$q_sensitivity<.05,na.rm=TRUE),
      direction_agreement=mean(paired$direction_agrees[good]),
      descriptive_coefficient_r=stats::cor(paired$coef_main[good],paired$coef_sensitivity[good]),
      scale_note=if(method=="CLR")"CLR and main effects have different scales" else "Threshold and model standardization settings follow their frozen analyses")
  }
  list(pairs=pub_bind(out),summary=pub_bind(summaries))
}

run_publication_tables <- function(b,run,config,export,attempt) {
  plans <- publication_table_plan(b,run,config)
  assert_true(length(plans)==config$publication$expected_workbooks,"Workbook plan must contain16 formal entries")
  index <- list()
  for (id in names(plans)) attempt(id,"table",function() {
    p <- plans[[id]]; if(is.function(p))p<-p()
    index[[id]] <<- pub_workbook(id,p$title,p$sheets,p$note,export,config,main=p$main)
    index[[id]]$content_id <<- submission_table_manifest()$content_id[match(id,submission_table_manifest()$table_id)]
  })
  safe_write_csv(submission_table_manifest(),file.path(export,"manifests/submission_table_identity.csv"))
  safe_write_csv(data.frame(content_id="medication_adjustment_sensitivity",previous_table_id="Supplementary_Table_S14",submission_table_id=NA_character_,status="INTERNAL_ARCHIVE_ONLY"),file.path(export,"manifests/retired_table_identity.csv"))
  safe_write_csv(submission_pathway_scope(run$root,b$candidates)$audit,file.path(export,"validation/S7_display_selection_audit.csv"))
  archived<-attr(plans,"internal_medication")
  for(n in names(archived$sheets))safe_write_csv(archived$sheets[[n]],file.path(export,"tables/archive/retired_medication_sensitivity",paste0(n,".csv")))
  # Detailed audit material is written separately from reader-facing workbooks.
  for(name in c("species","pooled","interactions","additive","additive_long"))
    safe_write_csv(b[[name]],file.path(export,"tables/archive/complete_statistical_columns",paste0(name,".csv")))
  for (module in intersect(c("metadata_audit","abundant_taxa","medication_sensitivity"),names(run$result_dirs))) {
    files<-list.files(file.path(run$result_dirs[[module]],"tables"),pattern="[.]csv$",full.names=TRUE,recursive=TRUE)
    for(f in files) {
      x<-pub_read_csv_snapshot(run,f)
      relative<-substring(f,nchar(file.path(run$result_dirs[[module]],"tables"))+2L)
      safe_write_csv(x,file.path(export,"tables/archive",module,relative))
    }
  }
  safe_write_csv(pub_bind(index),file.path(export,"table_index.csv"))
  pub_finalize_workbook_integrity(export)
  submission_validate_table_identity(pub_bind(index))
  invisible(index)
}
