#!/usr/bin/env Rscript
# Static/code/configuration checks only; never starts a model or opens a plot device.
local({
  if(.Platform$OS.type=="windows"&&!isTRUE(l10n_info()[["UTF-8"]]))Sys.setlocale("LC_CTYPE","English_United States.utf8")
  entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE));stopifnot(length(entry)==1L)
  root<-normalizePath(dirname(dirname(entry)),winslash="/",mustWork=TRUE)
  source(file.path(root,"R/functions/project_paths.R"),local=environment());activate_existing_project_library(root)
  source(file.path(root,"R/functions/publication_bootstrap.R"),local=environment());load_publication_support(root,environment())
  rows<-list();check<-function(name,action) {
    error<-"";pass<-tryCatch({action();TRUE},error=function(e){error<<-conditionMessage(e);FALSE})
    rows[[length(rows)+1L]]<<-data.frame(check=name,pass=pass,detail=error);message(if(pass)"PASS "else"FAIL ",name,if(nzchar(error))paste0(": ",error))
  }
  files<-unique(c(list.files(file.path(root,"R"),pattern="[.]R$",recursive=TRUE,full.names=TRUE),
    list.files(file.path(root,"tests"),pattern="[.]R$",recursive=TRUE,full.names=TRUE),list.files(root,pattern="[.]R$",full.names=TRUE)))
  check("All R syntax",function(){for(f in files)parse(f,encoding="UTF-8")})
  check("No legacy statistical input dependencies",function(){source(file.path(root,"tests/validate_no_legacy_dependency.R"),local=environment());validate_no_legacy_dependency(root)})
  check("Anonymous/network connections versus legacy filesystem guard",function(){
    source(file.path(root,"tests/test_release_access_guard.R"),local=environment())
    test_release_access_guard(root)
  })
  registry<-utils::read.csv(file.path(root,"config/pipeline_modules.csv"),stringsAsFactors=FALSE)
  check("LOCO registration and function dependency identities",function(){
    assert_true(!anyDuplicated(registry$module)&&registry$module[2]=="species_loco","LOCO must follow primary")
    for(f in registry$file)assert_file(file.path(root,"R/pipeline",f))
    source(file.path(root,"R/functions/analysis_bootstrap.R"),local=environment());load_analysis_support(root,environment())
    source(file.path(root,"R/pipeline/02a_species_cross_cohort_robustness.R"),local=environment())
    dependencies<-module_function_dependencies("species_loco","run_species_cross_cohort_robustness")
    assert_true(all(vapply(dependencies,exists,logical(1),envir=environment(),inherits=TRUE,mode="function")),"Missing LOCO checkpoint dependency")
    assert_true(!file.exists(file.path(root,"config/species_loco_candidates.csv")),"Candidate-only LOCO gate must not be present")
  })
  config<-read_config(root);config$project_root<-root;config$publication<-yaml::read_yaml(file.path(root,"config/publication_panel_config.yml"))
  manifest<-utils::read.csv(file.path(root,config$publication$panel_registry),stringsAsFactors=FALSE)
  check("Final publication manifest",function(){
    pub_validate_panel_counts(manifest,config)
    required_supplementary<-c("SFig1C_BMI25_cohort_alpha_combined","SFig3_taxon_pathway_assignment_Sankey",
      "SFig4A_sensitivity_BMI23_species","SFig4B_sensitivity_CLR_species","SFig4C_sensitivity_CLR_pathways")
    assert_true(all(required_supplementary%in%manifest$figure_id),"Missing retained supplementary panel")
    assert_true(!any(grepl("^SFig5|diversity_diagnostics",manifest$figure_id)),"Obsolete supplementary panel remains in the final manifest")
    assert_true(nrow(manifest)==40L&&sum(manifest$manuscript_figure=="Supplementary Figure 4")==3L,"Expected 40 final panels and three sensitivity scatters")
  })
  check("Migrated style function loading and latest heatmap geometry",function(){
    e<-final_style_context(list(),list(root=root),config,tempdir())
    required<-c("plot_alpha_diversity","plot_species_heatmap","plot_upset","plot_species_volcano","r3_plot_pathway_forest","r3_plot_pathway_sankey")
    assert_true(all(vapply(required,exists,logical(1),envir=e,inherits=TRUE,mode="function")),"Missing migrated visual dependency")
    x<-utils::read.csv(file.path(root,"config/refinement_round3_panels.csv"));z<-x[x$figure_id=="Fig6B_REML_pathway_effect_heatmap",]
    assert_true(z$width_mm==140&&z$height_mm==180,"Author's Fig6B geometry not applied")
  })
  check("Original and new statistical implementation hashes",function(){
    contract<-utils::read.csv(file.path(root,"config/statistical_implementation_sha256.csv"))
    assert_true(all(unname(vapply(file.path(root,contract$file),sha256_file,character(1)))==contract$sha256),"A statistical implementation contract changed")
  })
  check("Frozen package and resource descriptors",function(){
    lock<-jsonlite::fromJSON(file.path(root,"renv.lock"),simplifyVector=FALSE)
    assert_true(lock$Packages$Maaslin2$Version=="1.22.0"&&lock$Packages$curatedMetagenomicData$Version=="3.16.1"&&lock$Packages$colorspace$Version=="2.1-1","Frozen software descriptor mismatch")
    resource<-utils::read.csv(file.path(root,"config/source_resources.csv"));assert_true(!anyDuplicated(resource$resource_id)&&all(nchar(resource$sha256)==64L),"Resource identity/hash descriptor invalid")
  })
  check("Release validation plan covers current modules",function(){
    plan<-yaml::read_yaml(file.path(root,"config/release_validation_plan.yml"));ids<-vapply(plan$modules,function(x)x$module_id,character(1))
    assert_true(!anyDuplicated(ids)&&all(c("inputs",registry$module)%in%ids),"Release plan missing a computation module")
    admission<-yaml::read_yaml(file.path(root,"config/module_admission.yml"));assert_true(!isTRUE(admission$production_ready),"Code-only build cannot claim production admission")
  })
  out<-file.path(Sys.getenv("T2DBMI_AUDIT_DIR",unset=file.path(dirname(root),paste0(basename(root),"_Audit"))),paste0("StaticCheck_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_",Sys.getpid()));dir.create(out,recursive=TRUE)
  report<-do.call(rbind,rows);safe_write_csv(report,file.path(out,"static_checks.csv"))
  safe_write_lines(c(if(all(report$pass))"PASS"else"FAIL","Scope: syntax, paths, module/manifests, code hashes and dependency descriptors only. No model, plot/export, result comparison or clean-run acceptance was performed."),file.path(out,"scope.txt"))
  message("Static report: ",out);quit(status=if(all(report$pass))0L else 1L)
})
