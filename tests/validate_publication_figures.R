# Synthetic selections and a tiny TIFF fixture, not a study-figure export.
local({
  entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))
  if(!grepl("^([A-Za-z]:|/)",entry))entry<-file.path(getwd(),entry)
  root<-dirname(dirname(entry))
  source(file.path(root,"R/functions/project_paths.R"),local=environment());activate_existing_project_library(root)
  source(file.path(root,"R/functions/publication_bootstrap.R"),local=environment());load_publication_support(root,environment())
  config<-read_config(root);config$project_root<-root
  config$project_root<-root
  config$publication<-yaml::read_yaml(file.path(root,"config/publication_panel_config.yml"))
  config$output$font_family<-config$publication$font_family
  manifest<-utils::read.csv(file.path(root,config$publication$panel_registry))
  pub_validate_panel_counts(manifest,config)
  x<-data.frame(feature=paste0("species.Taxon.",sprintf("%02d",1:16)),
    coef=c(1:7,-(1:7),1,-1),qval=c(rep(.01,14),.051,.9))
  before<-serialize(x,NULL)
  a<-pub_select_volcano_labels(x)
  b<-pub_select_volcano_labels(x[nrow(x):1,,drop=FALSE])
  assert_true(sum(a$label_selected)==10L && sum(a$coef[a$label_selected]>0)==5L &&
    sum(a$coef[a$label_selected]<0)==5L,"Volcano direction caps failed")
  assert_true(identical(a$feature[order(a$label_rank,na.last=NA)],b$feature[order(b$label_rank,na.last=NA)]),
    "Label selection depends on incoming row order")
  assert_true(identical(serialize(x,NULL),before),"Label selection mutated source input")
  p<-pub_volcano(a,config,a[a$label_selected,])
  pub_validate_panel(p,a,"synthetic volcano",110,100)
  spec<-pub_comparisons()
  species<-pub_bind(lapply(1:5,function(i)data.frame(feature=paste0("species.Taxon.",1:200),
    comparison=spec$comparison[i],coef=seq(-1,1,length.out=200),
    qval=ifelse(1:200<=c(69,32,10,4,5)[i],.01,.5))))
  additive<-data.frame(feature=paste0("species.Taxon.",1:200),
    q_within_exposure_BH_T2D=ifelse(1:200<=62,.01,.5),
    q_within_exposure_BH_BMI_per5=ifelse(1:200%in%c(1:13,63:76),.01,.5))
  interaction<-pub_bind(lapply(c("BMI25","continuous_BMI","BMI23"),function(bmi)
    pub_bind(lapply(c("species","pathways"),function(domain)data.frame(BMI_definition=bmi,domain=domain,
      feature=seq_len(if(domain=="species")200L else 5945L),formal_q=.8)))))
  reml<-pub_bind(lapply(1:5,function(i){
    sig<-switch(as.character(i),`1`=1:13,`2`=c(1:4,14:142),`3`=integer(),`4`=integer(),`5`=141:160)
    positive<-switch(as.character(i),`1`=1:9,`2`=c(1:4,14:123),integer())
    data.frame(pathway=paste0("PWY",sprintf("%03d",1:200),":Synthetic description|g__Taxon.s__Taxon_example"),
      comparison=spec$comparison[i],padj=ifelse(1:200%in%sig,.01,.5),p_value=seq(.00001,.001,length.out=200),
      combined_log2FC=ifelse(1:200%in%positive,1,-1))
  }))
  fixture<-list(spec=spec,species=species,reml=reml,additive=additive,interactions=interaction)
  upset<-pub_species_upset(fixture,config)
  assert_true(sum(upset$intersections$n)==nrow(upset$data),"UpSet union does not reconcile")
  workflow<-pub_workflow_design(fixture,config)
  summary<-pub_graphical_summary(fixture,config)
  venn<-pub_overlap_venn(c(9,4,129),c("Overall","NT2D vs NC"),config)
  for(z in list(upset,workflow,summary,venn))pub_validate_panel(z$plot,z$data,"synthetic diagram",180,120)
  small<-pub_reml_heatmap_data(fixture,root)
  full<-pub_reml_heatmap_data(fixture,root,full=TRUE)
  assert_true(length(unique(small$pathway))<=35L && all(reml$pathway[reml$comparison==spec$comparison[1]&reml$padj<.05]%in%small$pathway),
    "Representative heatmap must include all overall significant features")
  assert_true(setequal(full$pathway,reml$pathway[reml$padj<.05]),"Full supplementary heatmap lost features")
  heat<-data.frame(comparison=c("A","B"),label=c("Taxon A","Taxon A"),coef=c(1,NA_real_),formal_q=c(.01,NA_real_))
  built<-ggplot2::ggplot_build(pub_heatmap(heat,config))
  assert_true(built$data[[1]]$fill[2]=="#E0E0E0" && identical(pub_stars(heat$formal_q),c("*","")),
    "Untested heatmap cells or stars changed")
  path<-tempfile(fileext=".tiff")
  grDevices::tiff(path,width=2,height=2,units="in",res=300,compression="lzw",bg="white",type="cairo")
  graphics::par(mar=c(0,0,0,0));graphics::plot.new();grDevices::dev.off()
  props<-pub_tiff_properties(path)
  assert_true(props$compression==5L && abs(props$x_dpi-300)<.01 && abs(props$y_dpi-300)<.01,
    "TIFF dpi/compression fixture failed")
  message("PASS: final-manifest panel registry, deterministic labels, UpSet/workflow/Venn/graphical-summary builds, representative/full heatmaps, NA stars and TIFF-format fixtures")
})
