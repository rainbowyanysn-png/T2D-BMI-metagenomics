run_species_cross_cohort_robustness <- function(data,run,config,expected) {
  require_project_packages(c("Maaslin2","tibble"));check_frozen_versions(config)
  assert_true(identical(config$models$primary_maaslin_standardize,TRUE) && config$filters$global_prevalence==.10 &&
    config$models$maaslin_normalization=="TSS" && config$models$maaslin_transform=="LOG","Frozen LOCO parameters changed")
  assert_true(nrow(data$metadata)==1186L && nrow(data$species)==200L && nrow(data$pathways)==5945L,"Frozen input sizes changed")
  assert_true(sum(data$metadata$disease=="T2D")==760L && sum(data$metadata$disease=="control")==426L,"Frozen disease counts changed")
  primary <- species_loco_primary_gate(run,expected)
  cohorts <- c("KarlssonFH_2013","MetaCardis_2020_a","QinJ_2012","SankaranarayananK_2015","YuJ_2015")
  assert_true(setequal(unique(as.character(data$metadata$study_name)),cohorts),"Frozen cohort IDs changed")
  significant <- species_loco_primary_significant(primary)
  assert_true(nrow(significant)==120L && length(unique(significant$feature))==79L,"Primary significant record/union checkpoint changed")
  started <- Sys.time()
  folder <- path_in_run(run,"tables","02a_species_LOCO",".keep");folder <- dirname(folder)
  save <- function(x,name) safe_write_csv(x,file.path(folder,name))
  save(significant,"LOCO_primary_significant_input.csv")
  save(data.frame(feature=rownames(data$species),order=seq_len(nrow(data$species))),"LOCO_global_feature_panel.csv")
  safe_write_lines(capture.output(sessionInfo()),file.path(folder,"sessionInfo.txt"))
  versions <- utils::read.csv(file.path(run$root,"config/package_versions.csv"))
  save(versions,"package_versions.csv")
  save(data.frame(comparison_id=names(primary),n_tested=vapply(primary,nrow,integer(1)),n_significant=vapply(primary,function(x)sum(x$qval<.05),integer(1)),status="PASS"),"LOCO_full_data_primary_gate.csv")
  all <- counts <- medications <- designs <- exclusions <- timings <- list(); k <- 0L
  specs <- species_loco_specs()
  set.seed(config$models$primary_seed)
  for(id in names(specs)) for(omitted in cohorts) {
    k <- k+1L; tag <- sprintf("L%02d",k); spec <- specs[[id]]; begin <- Sys.time()
    message("LOCO ",tag," ",id," omit ",omitted)
    result <- tryCatch({
      a <- species_loco_metadata(data,spec,omitted); m <- a$metadata
      x <- data$species[,m$sample_id,drop=FALSE]
      assert_true(identical(colnames(x),m$sample_id) && identical(rownames(x),rownames(data$species)),"LOCO sample/feature order changed")
      assert_true(nrow(x)==200L && length(unique(m$study_name))==4L,"LOCO must remove exactly one frozen cohort")
      expected_ids <- data$metadata$sample_id[as.character(data$metadata$study_name)!=omitted & as.character(data$metadata[[spec$column]]) %in% c(spec$reference,spec$comparison)]
      assert_true(identical(m$sample_id,expected_ids),"Unexpected LOCO sample exclusion/order")
      prevalence <- rowMeans(x>0); variable <- apply(x,1,stats::sd)>0
      # Maaslin2 1.22.0 uses STRICT >, not >=, for min_prevalence.
      eligible <- rowSums(x>0)>ncol(x)*.10 & variable
      reason <- ifelse(rowSums(x>0)<=ncol(x)*.10,"prevalence",ifelse(!variable,"zero_variance",""))
      model_dir <- dirname(path_in_run(run,"models",tag,".keep"))
      design_audit <- data.frame(model_id=tag,comparison_id=id,omitted_study=omitted,n_samples=nrow(m),n_reference=sum(m[[spec$variable]]==spec$reference),
        n_comparison=sum(m[[spec$variable]]==spec$comparison),n_cohorts=length(unique(m$study_name)),omitted_present=sum(m$study_name==omitted),
        design_rank=qr(a$design)$rank,design_columns=ncol(a$design),full_rank=qr(a$design)$rank==ncol(a$design),
        reference_group=spec$reference,comparison_group=spec$comparison,model_formula=paste("abundance ~",paste(spec$fixed,collapse=" + ")),
        references=paste(a$references,collapse=";"),medication_covariate=a$medication,medication_reference=a$medication_reference,
        medication_levels=if(a$medication)nlevels(m$antidiabetes_group) else NA_integer_,normalization="TSS",transform="LOG",standardize=TRUE,
        min_prevalence=.10,prevalence_comparator=">",seed=config$models$primary_seed,Maaslin2_version=as.character(utils::packageVersion("Maaslin2")),
        n_input_features=200L,n_eligible=sum(eligible),stringsAsFactors=FALSE)
      safe_write_csv(design_audit,file.path(model_dir,"preflight.csv"))
      fit <- audited_maaslin2(input_data=x,input_metadata=tibble::column_to_rownames(m,"sample_id"),output=model_dir,
        fixed_effects=spec$fixed,reference=a$references,normalization="TSS",transform="LOG",standardize=TRUE,min_prevalence=.10,
        plot_heatmap=FALSE,plot_scatter=FALSE)
      t <- fit$results[fit$results$metadata==spec$variable,,drop=FALSE]
      assert_true(!anyDuplicated(t$feature) && setequal(t$feature,rownames(x)[eligible]),"LOCO target panel differs from exact MaAsLin2 prevalence/variance rule")
      assert_true(all(t$value==spec$comparison) && all(t$N==nrow(m)),"LOCO target direction or sample count changed")
      assert_true(all(is.finite(t$coef)) && all(is.finite(t$stderr)) && all(is.finite(t$pval)) && all(is.finite(t$qval)),"Nonfinite LOCO target estimate")
      q_all <- stats::p.adjust(fit$results$pval,method="BH")
      assert_true(isTRUE(all.equal(fit$results$qval,q_all,tolerance=1e-14,check.attributes=FALSE)),"MaAsLin2 full-model BH q audit failed")
      t$qval_target_bh <- stats::p.adjust(t$pval,"BH")
      z <- data.frame(feature=rownames(x),tested_after_prevalence=eligible,prevalence=prevalence,exclusion_reason=reason,stringsAsFactors=FALSE)
      j <- match(z$feature,t$feature)
      for(n in c("metadata","value","coef","stderr","pval","qval","qval_target_bh","N","N.not.zero")) z[[n]] <- t[[n]][j]
      z$metadata <- spec$variable;z$value <- spec$comparison
      z$comparison_id <- id;z$omitted_study <- omitted;z$model_id <- tag
      z$ci_lower <- z$coef-1.96*z$stderr;z$ci_upper <- z$coef+1.96*z$stderr
      z$reference_group <- spec$reference;z$comparison_group <- spec$comparison
      z$model_formula <- design_audit$model_formula;z$medication_covariate <- a$medication;z$medication_reference <- a$medication_reference
      z$direction <- ifelse(is.na(z$coef),NA_character_,ifelse(z$coef>0,"Positive",ifelse(z$coef<0,"Negative","Zero")))
      z$q_family <- "MaAsLin2 full model qval";z$q_family_size <- sum(!is.na(fit$results$pval));z$n_target_tested <- nrow(t)
      z$normalization <- "TSS";z$transform <- "LOG";z$standardize <- TRUE;z$min_prevalence <- .10;z$Maaslin2_version <- design_audit$Maaslin2_version
      counts[[k]] <- data.frame(model_id=tag,comparison_id=id,omitted_study=omitted,study_name=as.character(m$study_name),sample_id=m$sample_id,group=as.character(m[[spec$variable]]))
      medications[[k]] <- if(a$medication) cbind(data.frame(model_id=tag,comparison_id=id,omitted_study=omitted),as.data.frame(table(group=m[[spec$variable]],antidiabetes_group=m$antidiabetes_group))) else NULL
      designs[[k]] <- design_audit;exclusions[[k]] <- z[!eligible,,drop=FALSE];all[[k]] <- z
      safe_write_csv(z,file.path(folder,paste0(tag,"_all_results.csv")))
      "PASS"
    },error=function(e) {
      safe_write_lines(c(id,omitted,conditionMessage(e)),file.path(folder,paste0(tag,"_FAIL.txt")))
      structure("FAIL",error=conditionMessage(e))
    })
    timings[[k]] <- data.frame(model_id=tag,comparison_id=id,omitted_study=omitted,started_utc=format(begin,tz="UTC",usetz=TRUE),
      finished_utc=format(Sys.time(),tz="UTC",usetz=TRUE),elapsed_seconds=as.numeric(difftime(Sys.time(),begin,units="secs")),status=as.character(result),error=if(is.null(attr(result,"error")))"" else attr(result,"error"))
    if(result!="PASS") {
      save(do.call(rbind,timings),"LOCO_run_summary_FAILED.csv")
      stop("LOCO model stopped without changing covariates: ",id," / ",omitted," / ",attr(result,"error"),call.=FALSE)
    }
  }
  results <- do.call(rbind,all);rownames(results)<-NULL
  sample_ids <- do.call(rbind,counts);rownames(sample_ids)<-NULL
  sample_counts <- stats::aggregate(sample_id~model_id+comparison_id+omitted_study+study_name+group,sample_ids,length)
  names(sample_counts)[names(sample_counts)=="sample_id"] <- "n"
  save(results,"LOCO_all_results.csv");save(do.call(rbind,designs),"LOCO_design_audit.csv")
  save(sample_ids,"LOCO_sample_ids.csv");save(sample_counts,"LOCO_sample_counts.csv");save(do.call(rbind,medications),"LOCO_medication_counts.csv")
  save(do.call(rbind,exclusions),"LOCO_prevalence_exclusions.csv");save(do.call(rbind,timings),"LOCO_run_summary.csv")
  summary <- species_loco_summarise_all(significant,results)
  features <- species_loco_feature_summary(summary)
  save(summary,"LOCO_primary_significant.csv")
  save(features,"LOCO_feature_summary.csv")
  species_loco_acceptance(summary,features)
  design <- do.call(rbind,designs);time <- do.call(rbind,timings)
  design$n_tested <- design$n_eligible
  design$n_significant <- vapply(design$model_id,function(id)sum(results$qval[results$model_id==id]<.05,na.rm=TRUE),integer(1))
  design$status <- time$status[match(design$model_id,time$model_id)]
  save(design,"LOCO_model_summary.csv")
  dictionary <- data.frame(field=names(results),definition=vapply(names(results),function(n)
    switch(n,qval="Formal MaAsLin2 full-model BH q; q<0.05 criterion",qval_target_bh="Audit only; target coefficient BH",tested_after_prevalence="Original MaAsLin2 strict >10% and nonzero variance; excluded estimates are NA",coef="Target comparison minus original reference coefficient",paste("Saved per-model",n)),character(1)))
  save(dictionary,"LOCO_data_dictionary.csv")
  save(data.frame(key=c("run_id","scope","n_loco_models","n_complete_input_features","n_all_result_rows","n_primary_significant_records","n_unique_significant_species","formal_q","elapsed_seconds"),
    value=c(run$id,"All five formal BMI25 species contrasts; medication adjusted T2D-internal model",25,200,nrow(results),nrow(summary),nrow(features),"MaAsLin2 full model qval",as.numeric(difftime(Sys.time(),started,units="secs")))),"run_manifest.csv")
  invisible(results)
}
