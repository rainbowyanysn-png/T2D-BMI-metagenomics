#!/usr/bin/env Rscript
# AUTHOR-INITIATED DIAGNOSTIC REFITS ONLY. Never called by run_all/export/validate.
local({
  if(.Platform$OS.type=="windows"&&!isTRUE(l10n_info()[["UTF-8"]]))Sys.setlocale("LC_CTYPE","English_United States.utf8")
  entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE));stopifnot(length(entry)==1L)
  root<-normalizePath(dirname(entry),winslash="/",mustWork=TRUE)
  source(file.path(root,"R/functions/project_paths.R"),local=environment());activate_existing_project_library(root)
  source(file.path(root,"R/functions/publication_bootstrap.R"),local=environment());load_publication_support(root,environment())
  source(file.path(root,"tests/release_validation_config.R"),local=environment())
  args<-commandArgs(TRUE);preflight_only<-"--preflight-only"%in%args;assert_true(sum(args=="--preflight-only")<=1L,"Duplicate --preflight-only");args<-args[args!="--preflight-only"];stopifnot(length(args)==4L,args[1]=="--run",args[3]=="--export",grepl("^Run_[0-9]{8}_[0-9]{6}$",args[2]),grepl("^PublicationExport_[0-9]{8}_[0-9]{6}([.][0-9]+)?_[0-9]+$",args[4]))
  run<-list(root=root,id=args[2],dir=file.path(root,"output",args[2]));run$result_dirs<-completed_result_directories(run)
  state<-latest_pipeline_state(run$dir)
  for(m in c("inputs","pathways","bmi23_pathways"))verify_module_checkpoint(run,checkpoint_record(run,state$completed[[m]]))
  config<-read_config(root);plan<-read_release_validation_plan(root)
  stopifnot(identical(config$models$pathway_meta_method,"REML"),requireNamespace("metafor",quietly=TRUE),requireNamespace("dplyr",quietly=TRUE))
  source(file.path(root,"tests/reml_diagnostic_environment.R"),local=environment())
  version_info<-diagnostic_metafor_version_check(root)
  inputs<-character();hashes<-character()
  read<-function(f){assert_file(f,"Diagnostic saved input");inputs<<-union(inputs,f);hashes[f]<<-sha256_file(f);utils::read.csv(f,check.names=FALSE,stringsAsFactors=FALSE)}
  source(file.path(root,"tests/reml_diagnostic_inputs.R"),local=environment())
  export<-file.path(root,"output",args[4])
  prepared<-diagnostic_prepare_inputs(root,run,export,read)
  forest_features<-prepared$forest_features
  # Both modes reach this point through the exact same complete real-input path.
  if(preflight_only){
    audit_root<-Sys.getenv("T2DBMI_AUDIT_DIR",unset=file.path(dirname(root),paste0(basename(root),"_Audit")))
    dir.create(audit_root,recursive=TRUE,showWarnings=FALSE)
    receipt<-file.path(audit_root,paste0("REML_Preflight_",invocation_stamp()));dir.create(receipt)
    safe_write_csv(prepared$audit,file.path(receipt,"input_preflight.csv"))
    actual<-vapply(inputs,sha256_file,character(1));unchanged<-actual==unname(hashes[inputs]);stopifnot(all(unchanged))
    safe_write_csv(data.frame(file=inputs,before_sha256=unname(hashes[inputs]),after_sha256=unname(actual),unchanged=unname(unchanged)),file.path(receipt,"read_only_input_hashes.csv"))
    safe_write_csv(version_info,file.path(receipt,"diagnostic_package_version_check.csv"))
    message("PASS complete diagnostic input preflight; comparisons: ",length(prepared$jobs),"; source files: ",length(inputs),". NO MODEL FITTED.\nPreflight report: ",receipt)
    quit(status=0L)
  }
  out<-file.path(root,"output",paste0("REML_Diagnostics_",invocation_stamp()));stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE)
  safe_write_csv(version_info,file.path(out,"diagnostic_package_version_check.csv"))
  safe_write_csv(prepared$audit,file.path(out,"input_preflight.csv"))
  message("AUTHOR DIAGNOSTIC REML REFITS; formal outputs remain read-only.\nDiagnostic directory: ",out)
  warning_rows<-all_results<-list();scalar<-function(x)if(length(x)==1L)as.numeric(x)else NA_real_
  for(job in prepared$jobs){
      m<-job$module;formal<-job$formal;comparison<-job$comparison;combined<-job$combined;paths<-job$paths
      records<-list()
      for(feature in paths){
        one<-if(m=="pathways")combined[combined$pathway==feature & !is.na(combined$lfcSE)&combined$lfcSE>0&is.finite(combined$log2FoldChange),,drop=FALSE]else combined[combined$pathway==feature & is.finite(combined$lfcSE)&combined$lfcSE>0&is.finite(combined$log2FoldChange),,drop=FALSE]
        if(nrow(one)<2L)next
        captured<-character();error<-""
        fit<-tryCatch(withCallingHandlers(metafor::rma(yi=log2FoldChange,sei=lfcSE,data=one,method=config$models$pathway_meta_method),warning=function(w){captured<<-c(captured,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e){error<<-conditionMessage(e);NULL})
        row<-data.frame(module=m,comparison=comparison,feature=feature,k=nrow(one),coef=if(is.null(fit))NA_real_ else scalar(fit$b),stderr=if(is.null(fit))NA_real_ else scalar(fit$se),pval=if(is.null(fit))NA_real_ else scalar(fit$pval),tau2=if(is.null(fit))NA_real_ else scalar(fit$tau2),I2=if(is.null(fit))NA_real_ else scalar(fit$I2),se_tau2=if(is.null(fit))NA_real_ else scalar(fit$se.tau2),warning_count=length(captured),error=error)
        j<-match(feature,formal$pathway);row$saved_coef<-formal$combined_log2FC[j];row$saved_stderr<-if("stderr"%in%names(formal))formal$stderr[j]else formal$se[j];row$saved_pval<-formal$p_value[j];row$saved_padj<-formal$padj[j];row$saved_tau2<-if("tau2"%in%names(formal))formal$tau2[j]else NA_real_
        row$coef_absolute_difference<-abs(row$coef-row$saved_coef);row$stderr_absolute_difference<-abs(row$stderr-row$saved_stderr);row$pval_absolute_difference<-abs(row$pval-row$saved_pval)
        row$finite_effect_and_se<-is.finite(row$coef)&&is.finite(row$stderr)&&row$stderr>0
        row$valid_nonnegative_tau2<-is.finite(row$tau2)&&row$tau2>=0
        row$within_original_absolute_tolerance<-row$coef_absolute_difference<=plan$absolute_tolerance & row$stderr_absolute_difference<=plan$absolute_tolerance & row$pval_absolute_difference<=plan$absolute_tolerance
        row$formal_significant<-row$saved_padj<.05
        row$entered_representative_forest<-m=="pathways"&&((comparison=="BMI25_Total_T2D_vs_HC"&&feature%in%forest_features$C1)||(comparison=="BMI25_NonOW_T2D_vs_NonOW_Ctrl"&&feature%in%forest_features$C2))
        row$discussion_status<-"MANUSCRIPT_NOT_SUPPLIED";row$diagnostic_fit_file<-NA_character_
        if(!is.null(fit)&&(length(captured)||!isTRUE(row$within_original_absolute_tolerance))){
          target<-diagnostic_fit_file(out,m,comparison,feature);dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
          assert_true(!file.exists(target),"Refusing to overwrite diagnostic fitted object (duplicate/collision)")
          saveRDS(fit,target);row$diagnostic_fit_file<-substring(target,nchar(out)+2L)
        }
        records[[feature]]<-row
        if(length(captured))for(w in captured)warning_rows[[length(warning_rows)+1L]]<-cbind(row,warning=w)
      }
      z<-dplyr::bind_rows(records);success<-!nzchar(z$error)
      z$diagnostic_q_BH<-NA_real_;z$diagnostic_family_matches_formal<-setequal(z$feature[success],formal$pathway)
      if(all(z$diagnostic_family_matches_formal))z$diagnostic_q_BH[success]<-p.adjust(z$pval[success],"BH")
      z$diagnostic_significance_matches_saved<-ifelse(is.na(z$diagnostic_q_BH),NA,(z$diagnostic_q_BH<.05)==(z$saved_padj<.05))
      safe_write_csv(z,file.path(out,paste0(comparison,"_diagnostic_refits.csv")));all_results[[paste(m,comparison)]]<-z
  }
  safe_write_csv(dplyr::bind_rows(all_results),file.path(out,"all_diagnostic_refits.csv"))
  safe_write_csv(if(length(warning_rows))dplyr::bind_rows(warning_rows)else data.frame(warning=character()),file.path(out,"feature_level_warnings.csv"))
  changed<-vapply(inputs,sha256_file,character(1))!=unname(hashes[inputs]);safe_write_csv(data.frame(file=inputs,before_sha256=unname(hashes[inputs]),after_sha256=vapply(inputs,sha256_file,character(1)),unchanged=!changed),file.path(out,"read_only_input_hashes.csv"));stopifnot(!any(changed))
  safe_write_lines(c("Diagnostic only. No formal result was replaced or release admitted.","Exact original rma REML formula/default optimizer and cohort scope, replayed from saved cohort effect/SE inputs. No DESeq2 or raw preprocessing rerun.","CSV serialization may differ slightly from original in-memory precision; warning counts and identity may not reproduce exactly. The original aggregate warning log remains unresolved until author review.","No alternative optimizer or profile likelihood was automatically substituted; finite SE/tau2 checks do not prove a global likelihood maximum. Review feature_level_warnings.csv and retained fitted objects before deciding whether a separately authorized profile/optimizer comparison is needed."),file.path(out,"README_DIAGNOSTIC_ONLY.txt"))
  safe_write_lines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"));message("Diagnostic complete; requires author review: ",out)
})
