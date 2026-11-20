# Migrated visual functions receive current-run objects, never external project files.
final_style_context <- function(b,run,config,export) {
  env<-new.env(parent=environment())
  env$r1_need<-assert_true;env$r3_need<-assert_true
  env$r2_columns<-function(d,columns)assert_columns(d,columns,"Migrated display data")
  env$r1_read<-function(root,paths,name) {
    if(name%in%c("Fig3C_intersection_counts.csv","Fig3C_set_sizes.csv")) {
      file<-file.path(export,"plot_data",name)
      assert_file(file);return(utils::read.csv(file,check.names=FALSE))
    }
    directory<-if(name=="abundant_taxa_emmeans_pairwise_all.csv")"04b_abundant_taxa"else if(name=="BMI25_alpha_LMM_summary.csv")"04_primary_diversity"else"04c_raw_Aitchison"
    pub_read(run,paste0("tables/",directory,"/",name))
  }
  files<-c("theme_round1","round1_annotations","plot_alpha_diversity","plot_ordination","plot_abundant_taxa","plot_cohort_permanova",
    "plot_composition","plot_upset","plot_lollipop","plot_species_heatmap","plot_species_volcano","plot_clinical_heatmap","pathway_pure_helpers","round3_renderers")
  for(f in files)source(file.path(run$root,"R/functions/refinement",paste0(f,".R")),local=env)
  env$style<-yaml::read_yaml(file.path(run$root,"config/refinement_style.yml"))
  env$pathway_style<-yaml::read_yaml(file.path(run$root,"config/refinement_pathway_style.yml"))
  env
}

final_cohort_alpha <- function(b,run,env) {
  d<-pub_read(run,"tables/04_primary_diversity/BMI25_alpha_per_sample.csv")
  d$group<-unname(group_display[as.character(d$analysis_group)])
  build<-function(metric) {
    z<-d;z$value<-z[[if(metric=="Shannon")"shannon"else"richness"]]
    panel<-data.frame(plot_family="alpha_cohort",metric=metric,title="Alpha Diversity by Cohort",width_mm=224)
    list(data=z,result=env$plot_alpha_diversity(z,panel,run$root,list(),env$style))
  }
  a<-build("Shannon");c<-build("Richness")
  assert_true(identical(a$data$sample_id,c$data$sample_id),"Combined alpha sample order changed")
  top<-a$result$plot+ggplot2::labs(title=NULL,x=NULL)+ggplot2::theme(axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank())
  bottom<-c$result$plot+ggplot2::labs(title=NULL,x="Cohort")+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=0,hjust=.5,vjust=1,margin=ggplot2::margin(t=2)))
  p<-patchwork::wrap_plots(top,bottom,ncol=1,guides="collect",heights=c(1,1)) & ggplot2::theme(legend.position="bottom",legend.direction="horizontal")
  p<-p+patchwork::plot_annotation(title="Alpha Diversity by Cohort",theme=ggplot2::theme(text=ggplot2::element_text(family=env$style$font$family),plot.title=ggplot2::element_text(size=env$style$font$title_size,face="bold",hjust=.5),plot.background=ggplot2::element_rect(fill="white",colour=NA)))
  a$data$metric_display<-"Shannon";c$data$metric_display<-"Richness"
  list(plot=p,data=pub_bind(list(a$data,c$data)),width=224,height=180)
}

final_panel_style <- function(id,title,data,plot,width,height,b,run,config,export,env) {
  registry<-utils::read.csv(file.path(run$root,"config/refinement_round1_panels.csv"),stringsAsFactors=FALSE)
  sfig<-utils::read.csv(file.path(run$root,"config/refinement_round1_sfig1_layout.csv"),stringsAsFactors=FALSE)
  r2<-utils::read.csv(file.path(run$root,"config/refinement_round2_panels.csv"),stringsAsFactors=FALSE)
  r3<-utils::read.csv(file.path(run$root,"config/refinement_round3_panels.csv"),stringsAsFactors=FALSE)
  result<-NULL;panel<-NULL;d<-data;caption_note<-NULL
  source_id<-switch(id,SFig1B_BMI25_richness="SFig1D_BMI25_richness",SFig1D_BMI25_ComBat_descriptive="SFig1F_BMI25_ComBat_descriptive",
    SFig1E_BMI25_ComBat_descriptive="SFig1G_BMI25_ComBat_descriptive",id)
  if(source_id%in%registry$figure_id) {
    panel<-registry[registry$figure_id==source_id,,drop=FALSE]
    if(id%in%sfig$output_filename) {
      row<-sfig[sfig$output_filename==id,,drop=FALSE];panel$title<-row$title;panel$width_mm<-row$width_mm;panel$height_mm<-row$height_mm
    }
    if(panel$plot_family%in%c("abundant_phylum","abundant_species")) {
      result<-env$plot_abundant_taxa(d,panel,run$root,list(),env$style)
      # Annotation rows include ALL 24/48 saved comparisons, not just visible brackets.
      original<-env$r1_read(run$root,list(),"abundant_taxa_emmeans_pairwise_all.csv")
      level<-if(panel$plot_family=="abundant_species")"species"else"phylum"
      original<-original[original$taxonomic_level==level,,drop=FALSE]
      for(n in names(original))assert_true(identical(original[[n]],result$statistics[[n]]),paste("Saved SFig2 annotation field changed:",n))
      shown<-result$statistics[result$statistics$display_selected,,drop=FALSE]
      layers<-Filter(function(x)inherits(x$geom,"GeomSegment")||inherits(x$geom,"GeomText"),result$plot$layers)
      assert_true(length(layers)==if(nrow(shown))4L else 0L,"Unexpected SFig2 annotation layers")
      for(layer in layers)for(n in c("feature","p_raw","p_adj_BH","display_x1","display_x2","display_y","display_label_y"))
        assert_true(identical(as.character(layer$data[[n]]),as.character(shown[[n]])),paste("Drawn SFig2 bracket differs from saved audit:",n))
      result$audits<-list(pairwise_annotations=result$statistics)
      safe_write_csv(result$statistics,file.path(export,"tables/archive",paste0(substr(id,1,6),"_pairwise_audit.csv")))
      safe_write_csv(result$statistics,file.path(export,"plot_data",paste0(substr(id,1,6),"_brackets.csv")))
      caption_note<-result$caption
    } else if(startsWith(panel$plot_family,"alpha")) {
      d$group<-if(panel$plot_family=="alpha_total")ifelse(d$disease=="control","Control","T2D")else unname(group_display[as.character(d$analysis_group)])
      result<-env$plot_alpha_diversity(d,panel,run$root,list(),env$style)
    } else if(startsWith(panel$plot_family,"combat")||startsWith(panel$plot_family,"ordination")) {
      # The original publication PCoA object stores phenotype but not disease/group metadata.
      # Recover those display labels by exact current-run sample ID, never by row position.
      k<-match(d$sample_id,b$md$sample_id)
      assert_true(!anyNA(k)&&!anyDuplicated(d$sample_id),"PCoA current-run metadata alignment failed")
      d$disease<-b$md$disease[k];d$analysis_group<-b$md$BMI25_group[k]
      d$phenotype<-if(panel$plot_family%in%c("combat_total","ordination_total"))ifelse(d$disease=="control","Control","T2D")else unname(group_display[as.character(d$analysis_group)])
      result<-env$plot_ordination(d,panel,run$root,list(),env$style)
      if(startsWith(panel$plot_family,"combat"))caption_note<-"Descriptive visualization only. Formal beta-diversity inference uses uncorrected CZM/CLR Aitchison distances."
    } else if(panel$plot_family=="cohort_permanova")result<-env$plot_cohort_permanova(d,panel,run$root,list(),env$style)
  } else if(id%in%r2$figure_id) {
    panel<-r2[r2$figure_id==id,,drop=FALSE]
    if(id=="Fig4F_Spearman_clinical_heatmap"&&!"feature_rank"%in%names(d))d$feature_rank<-b$clinical_selected$feature_rank[match(d$feature,b$clinical_selected$feature)]
    if(id=="Fig3C_BMI25_species_UpSet")for(n in b$spec$comparison)d[[n]]<-as.logical(d[[n]])
    fun<-switch(panel$plot_family,composition_phylum=env$plot_composition,composition_species=env$plot_composition,upset=env$plot_upset,
      lollipop=env$plot_lollipop,species_heatmap=env$plot_species_heatmap,volcano=env$plot_species_volcano,clinical_heatmap=env$plot_clinical_heatmap)
    result<-fun(d,panel,run$root,list(),env$style)
  } else if(id%in%r3$figure_id) {
    panel<-r3[r3$figure_id==id,,drop=FALSE]
    if(id%in%c("Fig5B_MaAsLin2_genus_pathway_features","Fig5C_MaAsLin2_species_pathway_features","SFig3_taxon_pathway_assignment_Sankey"))d<-b$candidates
    if(grepl("Fig6C[12]_",id)) {
      bundle<-final_pathway_forest(b,panel$panel);d<-bundle$values
      panel$title<-bundle$comparison
      safe_write_csv(bundle$ranking,file.path(export,"plot_data",paste0(id,"_candidate_ranking.csv")))
      safe_write_csv(bundle$pooled,file.path(export,"plot_data",paste0(id,"_pooled_source_row.csv")))
      safe_write_csv(bundle$cohorts,file.path(export,"plot_data",paste0(id,"_cohort_source_rows.csv")))
    }
    if(panel$plot_family=="pathway_overlap"&&!"set"%in%names(d))d$set<-c("Overall only","Shared","NT2D vs NC only")
    if(panel$plot_family=="pathway_categories"&&!"functional_category"%in%names(d))d$functional_category<-d$category
    result<-env$r3_build_panel(d,panel,env$pathway_style)
    if(panel$plot_family%in%c("pathway_taxonomy_bar","pathway_bubble","pathway_sankey")) {
      assert_true(identical(result$source_features,env$r3_expected_candidate_data(d)$feature),"Full pathway display family changed")
      caption_note<-"All 64 species-stratified pathway contributions meeting the prespecified visualization criteria are displayed (six species, four genera and 42 pathway IDs); no top-N truncation. Assignments represent inferred functional potential, not metabolic flux or a causal mechanism."
    }
  }
  if(!is.null(result)) {
    plot<-result$plot;data<-result$display_data;width<-panel$width_mm;height<-panel$height_mm
    if(length(result$audits))for(audit_name in names(result$audits)) {
      n<-sub("[.]csv$","",audit_name)
      safe_write_csv(result$audits[[audit_name]],file.path(export,"plot_data",paste0(id,"_display_",n,".csv")))
    }
  }
  if(startsWith(id,"SFig4") && grepl("_sensitivity_",id,fixed=TRUE)) {
    # Preserve the sensitivity-specific larger typography from the shared builder.
    plot<-plot+ggplot2::theme(legend.background=ggplot2::element_blank(),
      legend.box.background=ggplot2::element_blank(),legend.key=ggplot2::element_blank(),axis.text.y=ggplot2::element_text(face="plain"))
  }
  list(plot=plot,data=data,width=width,height=height,note=caption_note,title=if(!is.null(panel))panel$title else NULL)
}

# Shared display-only builder for full export and single-part resizing.
# Uses the existing saved display classes and facet text; no summaries/tests are computed.
# Display-only: move existing four saved annotation lines into facet strips.
final_sensitivity_strip_labels <- function(annotation) {
  assert_columns(annotation,c("comparison","label"),"Saved sensitivity annotations")
  assert_true(!anyDuplicated(annotation$comparison),"Duplicate annotation facet")
  lines<-strsplit(as.character(annotation$label),"\n",fixed=FALSE)
  assert_true(all(lengths(lines)==4L),"Expected four existing saved descriptive annotation lines")
  labels<-vapply(seq_len(nrow(annotation)),function(i)paste0(annotation$comparison[i],"\n",lines[[i]][1],"; ",lines[[i]][2],"\n",lines[[i]][3],"; ",lines[[i]][4]),character(1))
  setNames(labels,as.character(annotation$comparison))
}

final_sensitivity_scatter_plot <- function(d,annotation,spec,config) {
  assert_columns(d,c("coef_main","coef_sensitivity","comparison","direction_class","shared_status","shared_fdr"),"Saved sensitivity plot-data")
  assert_columns(annotation,c("comparison","label"),"Saved sensitivity facet annotations")
  assert_true(!anyNA(d$shared_fdr)&&all(is.finite(d$coef_main))&&all(is.finite(d$coef_sensitivity)),"Saved sensitivity display has missing coefficients/classes")
  assert_true(all(as.character(d$direction_class)%in%c("Concordant","Discordant"))&&all(as.character(d$shared_status)%in%c("Not shared FDR-significant","Shared FDR-significant")),"Saved sensitivity display classes changed")
  assert_true(!anyDuplicated(annotation$comparison)&&setequal(annotation$comparison,unique(d$comparison)),"Saved facet annotation coverage differs")
  d$direction_class<-factor(as.character(d$direction_class),levels=c("Concordant","Discordant"))
  d$shared_status<-factor(as.character(d$shared_status),levels=c("Not shared FDR-significant","Shared FDR-significant"))
  strip_annotations<-identical(spec$id,"SFig4C_sensitivity_CLR_pathways")
  strip_labels<-if(strip_annotations)final_sensitivity_strip_labels(annotation)else NULL
  point_size<-if(spec$domain=="pathways").45 else .65
  direction_colours<-c(Concordant="#2166AC",Discordant="#F7AA58")
  shared_alpha<-c("Not shared FDR-significant"=.18,"Shared FDR-significant"=.95)
  typography<-config$publication$sensitivity_scatter
  required<-c("base_size_pt","title_size_pt","axis_title_size_pt","axis_text_size_pt","strip_text_size_pt","legend_title_size_pt","legend_text_size_pt","annotation_size_mm")
  values<-unlist(typography[required])
  assert_true(length(values)==length(required)&&all(is.finite(values))&&all(values>0),"Missing/invalid SFig4 typography settings")
  display_config<-config;display_config$publication$font_size<-typography$base_size_pt
  theme<-publication_panel_theme(display_config)+ggplot2::theme(
    text=ggplot2::element_text(family=config$output$font_family,size=typography$base_size_pt),
    plot.title=ggplot2::element_text(size=typography$title_size_pt,face="bold"),
    axis.title=ggplot2::element_text(size=typography$axis_title_size_pt,face="bold"),
    axis.text=ggplot2::element_text(size=typography$axis_text_size_pt,colour="black"),
    strip.text=ggplot2::element_text(size=typography$strip_text_size_pt),
    legend.title=ggplot2::element_text(size=typography$legend_title_size_pt),
    legend.text=ggplot2::element_text(size=typography$legend_text_size_pt))
      plot<-ggplot2::ggplot(d,ggplot2::aes(coef_main,coef_sensitivity,colour=direction_class,alpha=shared_status))+
        ggplot2::geom_point(data=d[!d$shared_fdr,,drop=FALSE],shape=16,size=point_size)+
        ggplot2::geom_point(data=d[d$shared_fdr,,drop=FALSE],shape=16,size=point_size)+
        ggplot2::geom_hline(yintercept=0,colour="grey75",linewidth=.25)+
        ggplot2::geom_vline(xintercept=0,colour="grey75",linewidth=.25)+
        ggplot2::facet_wrap(~comparison,scales="free",labeller=if(strip_annotations)ggplot2::as_labeller(strip_labels)else ggplot2::label_value)+
        ggplot2::scale_colour_manual(values=direction_colours,drop=FALSE)+
        ggplot2::scale_alpha_manual(values=shared_alpha,drop=FALSE)+
        ggplot2::labs(title=spec$title,x=spec$xlab,y=spec$ylab,colour="Coefficient direction",alpha="Shared FDR significance")+
        theme+ggplot2::theme(legend.background=ggplot2::element_blank(),legend.box.background=ggplot2::element_blank(),
          legend.key=ggplot2::element_blank(),axis.text.y=ggplot2::element_text(face="plain"))
  if(strip_annotations) {
    # Facet strips are outside the data panel, so text cannot obscure points.
    plot<-plot+ggplot2::theme(strip.placement="outside",strip.background=ggplot2::element_rect(fill="white",colour=NA),
      strip.text=ggplot2::element_text(size=typography$annotation_size_mm*72.27/25.4,lineheight=1.05,
        margin=ggplot2::margin(t=2,r=2,b=3,l=2),hjust=.5))
  } else plot<-plot+ggplot2::geom_text(data=annotation,ggplot2::aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,
    family=config$output$font_family,hjust=-.05,vjust=1.10,size=typography$annotation_size_mm,lineheight=.95)
  plot
}
