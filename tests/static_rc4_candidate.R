#!/usr/bin/env Rscript
# RC4 STATIC ACCEPTANCE ONLY: no models, downloads, statistical export or plot devices.
local({
  if(.Platform$OS.type=="windows"&&!isTRUE(l10n_info()[["UTF-8"]]))Sys.setlocale("LC_CTYPE","English_United States.utf8")
  entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE));stopifnot(length(entry)==1L)
  root<-normalizePath(dirname(dirname(entry)),winslash="/",mustWork=TRUE)
  source(file.path(root,"R/functions/project_paths.R"),local=environment());activate_existing_project_library(root)
  source(file.path(root,"R/functions/publication_bootstrap.R"),local=environment());load_publication_support(root,environment())
  source(file.path(root,"R/functions/public_release_files.R"),local=environment())
  results<-list();check<-function(n,expr){detail<-"";ok<-tryCatch({force(expr);TRUE},error=function(e){detail<<-conditionMessage(e);FALSE});results[[length(results)+1L]]<<-data.frame(check=n,pass=ok,detail=detail);message(if(ok)"PASS "else"FAIL ",n,if(nzchar(detail))paste(":",detail))}
  text<-function(n)paste(readLines(file.path(root,n),encoding="UTF-8",warn=FALSE),collapse="\n")
  check("All R syntax",{files<-unique(c(list.files(file.path(root,"R"),pattern="[.]R$",recursive=TRUE,full.names=TRUE),list.files(file.path(root,"tests"),pattern="[.]R$",recursive=TRUE,full.names=TRUE),list.files(root,pattern="[.]R$",full.names=TRUE)));for(f in files)parse(f,encoding="UTF-8")})
  check("31 frozen statistical hashes",{z<-utils::read.csv(file.path(root,"config/statistical_implementation_sha256.csv"));stopifnot(nrow(z)==31,all(vapply(file.path(root,z$file),sha256_file,character(1))==z$sha256))})
  registry<-utils::read.csv(file.path(root,"config/pipeline_modules.csv"),stringsAsFactors=FALSE)
  check("Inputs plus 22 registered model/preparation modules and dependencies",{
    stopifnot(nrow(registry)==22,!anyDuplicated(registry$module),!any(registry$module%in%c("extreme","power","core_score","fb_ratio","pilot")))
    e<-new.env(parent=environment());source(file.path(root,"R/functions/analysis_bootstrap.R"),local=e);e$load_analysis_support(root,e)
    for(i in seq_len(nrow(registry))){f<-file.path(root,"R/pipeline",registry$file[i]);stopifnot(file.exists(f));source(f,local=e);stopifnot(exists(registry$function_name[i],envir=e,mode="function"));for(n in module_function_dependencies(registry$module[i],registry$function_name[i]))stopifnot(exists(n,envir=e,inherits=TRUE,mode="function"))}
    stopifnot(file.exists(file.path(root,"R/pipeline/01_import_harmonize_qc.R")))
  })
  check("Pending RC4 admission and frozen environment",{
    ad<-yaml::read_yaml(file.path(root,"config/module_admission.yml"));stopifnot(!isTRUE(ad$production_ready),!isTRUE(ad$clean_run_pass),isTRUE(ad$clean_rc4_run_pending))
    lock<-jsonlite::fromJSON(file.path(root,"renv.lock"),simplifyVector=FALSE);stopifnot(lock$R$Version=="4.5.1",lock$Packages$Maaslin2$Version=="1.22.0",lock$Packages$curatedMetagenomicData$Version=="3.16.1")
  })
  check("RC4 creation and old-Run rejection",{
    sandbox<-tempfile("rc4_identity_");dir.create(sandbox);on.exit(unlink(sandbox,recursive=TRUE),add=TRUE)
    run<-run_context(sandbox);assert_rc4_run_identity(run);stopifnot(!file.exists(file.path(run$dir,"RC3_RUN_IDENTITY.txt")))
    bad<-list(dir=file.path(sandbox,"foreign"));dir.create(bad$dir);writeLines("RC3_CURRENT_PROJECT_RUN",file.path(bad$dir,"RC3_RUN_IDENTITY.txt"))
    stopifnot(inherits(tryCatch({assert_rc4_run_identity(bad);NULL},error=function(e)e),"error"))
    stopifnot(grepl("assert_rc4_run_identity(run)",text("R/functions/pipeline_runner.R"),fixed=TRUE),grepl("assert_rc4_run_identity(run)",text("R/functions/checkpoints.R"),fixed=TRUE))
  })
  check("No legacy runtime dependency or hardcoded source path",{source(file.path(root,"tests/validate_no_legacy_dependency.R"),local=environment());validate_no_legacy_dependency(root)})
  check("Legacy path guard; downloads and anonymous connections remain allowed",{
    source(file.path(root,"tests/test_release_access_guard.R"),local=environment());test_release_access_guard(root)
    e<-new.env(parent=environment());source(file.path(root,"tests/release_access_guard.R"),local=e);e$rc3_install_guard(root)
    tryCatch({
      g<-get("v3_guard",.GlobalEnv);g("");g("https://experimenthub.bioconductor.org");g(file.path(root,"config/analysis_config.yml"))
      for(n in c("T2D_BMI_RC3_FinalIntegration","T2D_BMI_RC3_CleanTest","T2D_BMI_RC3_IntegrationPreview","T2D_BMI_RC3_FigureRefinement","T2D_BMI_RC3_LOCO_Trial","T2D_BMI_Publication_Package_20261113_RC2","1015返修分析后文件整理V3","1115返修分析后文件整理V3_MedicationRestored")){
        err<-tryCatch({g(file.path(dirname(root),n,"output/result.csv"));NULL},error=function(e)e);stopifnot(inherits(err,"error"))
      }
    },finally={for(n in c("file","gzfile","bzfile","xzfile","unz","readRDS","load","source","sys.source","readLines","scan"))suppressMessages(untrace(n,where=baseenv()));for(n in c("v3_guard","v3_file_access_report"))if(exists(n,.GlobalEnv,inherits=FALSE))rm(list=n,envir=.GlobalEnv)})
  })
  config<-read_config(root);config$project_root<-root;config$publication<-yaml::read_yaml(file.path(root,"config/publication_panel_config.yml"))
  check("40 panels, 16 workbooks; latest drawing settings",{
    m<-utils::read.csv(file.path(root,"config/publication_panels_manifest.csv"));pub_validate_panel_counts(m,config);stopifnot(nrow(m)==40,config$publication$expected_workbooks==16,!any(grepl("^SFig[56]",m$figure_id)))
    r1<-utils::read.csv(file.path(root,"config/refinement_round1_panels.csv"));z<-r1[startsWith(r1$figure_id,"SFig2"),];stopifnot(all(z$width_mm==180),identical(z$height_mm,c(47.5,90)))
    r3<-utils::read.csv(file.path(root,"config/refinement_round3_panels.csv"));z<-r3[r3$figure_id=="Fig6B_REML_pathway_effect_heatmap",];stopifnot(z$width_mm==140,z$height_mm==180)
    z<-r3[r3$figure_id=="SFig3_taxon_pathway_assignment_Sankey",];stopifnot(z$width_mm==240,z$height_mm==180,config$publication$tiff_dpi_by_figure$SFig3_taxon_pathway_assignment_Sankey==600,config$publication$sensitivity_scatter$width_mm==234,config$publication$sensitivity_scatter$height_mm==115)
    a<-text("R/functions/publication_refinement_adapter.R");stopifnot(grepl("result<-env$plot_abundant_taxa(",a,fixed=TRUE),!grepl('!startsWith(id,"SFig2")',a,fixed=TRUE),grepl('guides="collect",heights=c(1,1)',a,fixed=TRUE))
    f<-text("R/functions/pathway_forest_selection.R");stopifnot(grepl("P162-PWY",f,fixed=TRUE),grepl("VALSYN-PWY",f,fixed=TRUE),grepl("nrow(one)==5L",f,fixed=TRUE))
    stopifnot(grepl("full 64",text("R/functions/refinement/round3_renderers.R"),fixed=TRUE)||grepl("All 64",text("R/functions/refinement/round3_renderers.R"),fixed=TRUE))
  })
  check("S3/S4/LOCO/new-S14 coverage and current-run source routing",{
    tables<-text("R/pipeline/12b_publication_tables.R");stopifnot(grepl("Supplementary_Table_S14",tables,fixed=TRUE),grepl("Database_completeness",tables,fixed=TRUE),grepl("Categorical antidiabetic treatment was included",tables,fixed=TRUE),!grepl("S15_REVIEW",tables,fixed=TRUE))
    f<-text("R/functions/final_integration_tables.R");for(n in c("LOCO_primary_significant","LOCO_feature_summary","LOCO_all_results","LOCO_model_summary","LOCO_exclusions","LOCO_data_dictionary"))stopifnot(grepl(n,f,fixed=TRUE))
    plan<-yaml::read_yaml(file.path(root,"config/release_validation_plan.yml"));stopifnot(all(c("inputs",registry$module)%in%vapply(plan$modules,function(z)z$module_id,character(1))))
    d<-pub_source_review(root,config,"confounder");stopifnot(nrow(d)==75,identical(pub_reader_table(d)$reason_not_adjusted,d$reason_not_adjusted))
  })
  check("Empty output/derived; no preview scripts or cached results",{
    for(n in c("output","data_derived")){files<-list.files(file.path(root,n),recursive=TRUE,all.files=TRUE,no..=TRUE);stopifnot(all(files==".gitkeep"))}
    stopifnot(!dir.exists(file.path(root,"audit")),!file.exists(file.path(root,"preview_publication.R")),!file.exists(file.path(root,"preview_s15.R")))
  })
  check("Public payload selected by shared policy and exact hashes",{
    x<-public_release_selection(root);saved<-utils::read.csv(file.path(root,"public_file_manifest.csv"),stringsAsFactors=FALSE);stopifnot(identical(x$file,saved$file),identical(x$sha256,saved$sha256),identical(x$bytes,as.numeric(saved$bytes)))
  })
  folder<-Sys.getenv("T2DBMI_AUDIT_DIR",unset=file.path(dirname(root),paste0(basename(root),"_Audit")))
  dir.create(folder,recursive=TRUE,showWarnings=FALSE);out<-file.path(folder,paste0("RC4_Static_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_",Sys.getpid(),".csv"))
  utils::write.csv(do.call(rbind,results),out,row.names=FALSE);message("Static report: ",out)
  if(!all(vapply(results,function(z)z$pass,logical(1))))quit(status=1L)
})
