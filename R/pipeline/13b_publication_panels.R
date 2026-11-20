run_publication_panels <- function(b,run,config,export,attempt) {
  refinement<-final_style_context(b,run,config,export)
  registry<-utils::read.csv(file.path(run$root,config$publication$panel_registry),stringsAsFactors=FALSE)
  index<-list();pending<-list();preflight<-list();theme<-publication_panel_theme(config);cols<-pub_effect_colours(config)
  style<-read_presentation_config(run$root);groups<-unlist(style$groups$fill)
  read<-function(p)pub_read(run,paste0("tables/",p))
  save<-function(id,fig,panel,title,plot,data,caption,rule,q="Not applicable",model="Descriptive",width=NULL,height=NULL,supp=FALSE){
    id<-switch(id,SFig1D_BMI25_richness="SFig1B_BMI25_richness",SFig1F_BMI25_ComBat_descriptive="SFig1D_BMI25_ComBat_descriptive",SFig1G_BMI25_ComBat_descriptive="SFig1E_BMI25_ComBat_descriptive",id)
    if(!id%in%registry$figure_id)return(invisible(NULL))
    registered<-registry[registry$figure_id==id,,drop=FALSE];fig<-registered$manuscript_figure;panel<-registered$panel
    migrated<-final_panel_style(id,title,data,plot,width,height,b,run,config,export,refinement)
    plot<-migrated$plot;data<-migrated$data;width<-migrated$width;height<-migrated$height
    if(!is.null(migrated$title))title<-migrated$title
    if(id%in%c("Fig5B_MaAsLin2_genus_pathway_features","Fig5C_MaAsLin2_species_pathway_features","SFig3_taxon_pathway_assignment_Sankey")) {
      caption<-"Species-stratified pathway features meeting the prespecified visualization criterion: full-model q<0.01 and absolute coefficient >1, with named taxonomic contributors; UNINTEGRATED/UNMAPPED excluded. All 64 eligible assignments are displayed; taxonomic assignments indicate inferred functional potential, not metabolic flux or a causal mechanism."
      rule<-"All eligible species-stratified assignments; full-model q<0.01 and absolute coefficient >1; ascending q, descending absolute coefficient, feature ID; no top-N truncation"
    }
    if(!is.null(migrated$note))caption<-paste(caption,migrated$note)
    if(id%in%c("Fig2C_BMI25_raw_Aitchison","Fig2D_BMI25_raw_Aitchison"))caption<-paste(caption,final_beta_caption(run))
    assert_true(is.null(pending[[id]]),paste("Duplicate panel ID:",id))
    width<-if(is.null(width))config$publication$half_width_mm else width
    height<-if(is.null(height))config$publication$default_height_mm else height
    title<-pub_caption_text(title);caption<-pub_caption_text(caption);rule<-pub_caption_text(rule)
    pub_language_check(c(title,caption,rule),id)
    pub_validate_panel(plot,data,id,width,height)
    pending[[id]]<<-list(plot=plot,id=id,figure=fig,panel=panel,title=title,data=data,caption=caption,
      rule=rule,export=export,config=config,q=q,model=model,width=width,height=height,supp=supp)
  }
  plan<-function(id,action){
    if(id%in%c("Fig2F","SFig6"))return(invisible(TRUE))
    ok<-attempt(id,"figure_preparation",action)
    preflight[[id]]<<-data.frame(task=id,preparation_complete=ok)
    invisible(ok)
  }
  basecaption<-"BMI25: NC, non-overweight control; OC, overweight control; NT2D, non-overweight T2D; OT2D, overweight T2D; BMI <25 versus >=25 kg/m²."

  plan("Fig1_workflow",function(){
    z<-pub_workflow_design(b,config)
    save("Fig1_workflow","Figure 1","","Study design and analysis workflow",z$plot,z$data,
      paste(basecaption,"Primary, additional and sensitivity analyses are parallel branches of the prepared data; no analysis is inferred from a preceding result branch."),
      "Six modules; sample counts from saved metadata",width=180,height=140)
  })

  alpha<-read("04_primary_diversity/BMI25_alpha_per_sample.csv")
  asum<-read("04_primary_diversity/BMI25_alpha_LMM_summary.csv")
  for(metric in c("shannon","richness"))for(scope in c("total","four_groups","cohort"))local({
    metric<-metric;scope<-scope
    id<-if(metric=="shannon" && scope!="cohort")paste0("Fig2",if(scope=="total")"A" else "B","_BMI25_Shannon") else
      paste0("SFig1",if(metric=="shannon")"B" else switch(scope,total="A",four_groups="D",cohort="C"),"_BMI25_",metric)
    plan(id,function(){
      d<-alpha;d$value<-d[[metric]];d$group<-if(scope=="total")factor(ifelse(d$disease=="control","Control","T2D"),levels=c("Control","T2D")) else
        factor(group_display[d$analysis_group],levels=names(groups))
      palette<-if(scope=="total")unlist(style$groups$total) else groups
      if(scope=="total")names(palette)[names(palette)=="Ctrl"]<-"Control"
      p<-ggplot2::ggplot(d,ggplot2::aes(group,value,fill=group))+
        ggplot2::geom_boxplot(width=.5,outlier.shape=NA,linewidth=.35,alpha=.7)+
        ggplot2::geom_point(position=ggplot2::position_jitter(width=.12,seed=config$output$label_seed),size=.35,alpha=.3)+
        ggplot2::scale_fill_manual(values=pub_palette(palette,d$group,paste(id,"groups")))+ggplot2::labs(x=NULL,y=if(metric=="shannon")"Shannon index" else "Observed species richness")+
        theme+ggplot2::theme(legend.position="none")
      if(scope=="cohort")p<-p+ggplot2::facet_wrap(~study_name,nrow=1)+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1))
      row<-asum[asum$Metric==if(metric=="shannon")"Shannon" else "Richness",]
      caption<-paste(basecaption,"LMM results. Overall T2D comparison: P",saved_p_label(row$Overall_T2D_vs_Control_p),
        "; four-group comparison: P",saved_p_label(row$Four_Group_p),"; ICC =",row$ICC,
        ". Boxplots show median and interquartile range; points show participants.")
      save(id,if(startsWith(id,"Fig2"))"Figure 2" else "Supplementary Figure 1",sub("^.*Fig[12]([A-D]).*$","\\1",id),
        paste("BMI25",metric,"diversity"),p,d,caption,"All current-run alpha-diversity observations; no model refitting",
        model="Frozen alpha LMM",width=if(scope=="cohort")180 else 90,height=if(scope=="cohort")75 else 85,supp=startsWith(id,"SFig"))
    })
  })
  for(scope in c("total","four_groups"))local({scope<-scope;panel<-if(scope=="total")"C" else "D"
    plan(paste0("Fig2",panel),function(){
      coords<-read("04c_raw_Aitchison/BMI25_raw_Aitchison_PCoA_coordinates.csv")
      md<-read("04c_raw_Aitchison/BMI25_raw_Aitchison_sample_manifest.csv")
      z<-pub_pcoa(coords,md,config,total=scope=="total")
      caption<-paste(basecaption,"Uncorrected CZM/CLR Aitchison distances;",nrow(coords),"samples. Dashed 95% t data ellipses describe within-group distributions.",
        "Four-group PERMDISP is significant; ordination separation cannot be attributed entirely to group centroid differences. Formal tests are provided separately.")
      save(paste0("Fig2",panel,"_BMI25_raw_Aitchison"),"Figure 2",panel,"Uncorrected Aitchison ordination",z$plot,z$data,
        caption,"Saved current-run PCoA coordinates; 95% t ellipses",model="Uncorrected Aitchison PCoA",width=110,height=85)
      geometry<-ggplot2::ggplot_build(z$plot)$data[[2]]
      safe_write_csv(geometry,file.path(export,"plot_data",paste0("Fig2",panel,"_ellipse_geometry.csv")))

    })
  })
  plan("Fig2E",function(){
    d<-read("04c_raw_Aitchison/BMI25_raw_Aitchison_PERMANOVA_by_cohort.csv")
    z<-pub_diagnostic_plot(d,config)
    save("Fig2E_BMI25_cohort_PERMANOVA","Figure 2","E","Cohort-specific PERMANOVA",z$plot,z$data,
      "Uncorrected Aitchison distances; four-group marginal tests within cohorts, 999 permutations. Model formulas and omitted nonvarying nuisance terms accompany plot-data.",
      "All five saved cohort tests",model="Raw Aitchison marginal PERMANOVA",width=150,height=75)

  })
  abundant<-read("04b_abundant_taxa/abundant_taxa_sample_level_data.csv")
  for(rank in c("phylum","species"))local({rank<-rank;panel<-if(rank=="phylum")"A" else "B"
    plan(paste0("Fig3",panel),function(){
      # Plot only the separately prepared and saved current-run display table.
      # 'Other' is the sum of actual omitted taxa, never a forced complement.
      d<-read(paste0("01c_publication_display/BMI25_",rank,"_composition_sample_level.csv"))
      assert_columns(d,c("sample_id","analysis_group","feature","raw_abundance","abundance_unit"),"Composition display data")
      stopifnot(all(is.finite(d$raw_abundance)),all(d$raw_abundance>=0),
        setequal(d$sample_id,b$md$sample_id),!anyDuplicated(d[,c("sample_id","feature")]))
      d$group<-factor(display_group(d$analysis_group),levels=names(groups))
      d$scope<-"BMI25 four groups"
      overall<-d;overall$group<-ifelse(b$md$disease[match(overall$sample_id,b$md$sample_id)]=="control","Control","T2D")
      overall$scope<-"Overall disease groups"
      m<-stats::aggregate(raw_abundance~scope+group+feature,pub_bind(list(overall,d)),mean)
      m$scope<-factor(m$scope,levels=c("Overall disease groups","BMI25 four groups"))
      m$group<-factor(m$group,levels=c("Control","T2D","NC","OC","NT2D","OT2D"))
      m$display_name<-if(rank=="species")pub_name(m$feature) else m$feature
      if(rank=="phylum"){
        palette<-unlist(config$publication$phylum_colours)
        palette<-c(palette[names(palette)!="Other"],`Other phyla`=palette[["Other"]])
      } else {
        selected<-read("04b_abundant_taxa/selected_top8_species.csv")
        selected<-selected[order(selected$display_order),]
        labels<-c(pub_name(selected$feature),"Other retained features")
        palette<-setNames(unlist(config$publication$species_palette),labels)
      }
      stopifnot(all(unique(m$display_name)%in%names(palette)))
      m$display_name<-factor(m$display_name,levels=names(palette))
      p<-ggplot2::ggplot(m,ggplot2::aes(group,raw_abundance,fill=display_name))+ggplot2::geom_col(width=.65,colour="white",linewidth=.15)+
        ggplot2::facet_grid(~scope,scales="free_x",space="free_x")+
        ggplot2::scale_fill_manual(values=pub_palette(palette,m$display_name,paste("Fig3",panel,"taxa")))+ggplot2::labs(x=NULL,y=if(rank=="species")"Mean relative abundance within the retained species panel" else "Mean relative abundance (%)",fill=if(rank=="species")"Species" else "Phylum")+theme
      if(rank=="species")p<-p+ggplot2::theme(legend.text=ggplot2::element_text(face="italic"))
      save(paste0("Fig3",panel,"_BMI25_",rank,"_composition"),"Figure 3",panel,paste("BMI25",rank,"composition"),p,m,
        paste(basecaption,if(rank=="species")"Top 8 species and remaining retained features; no reclosure to 100%." else "Top 4 phyla and the actual sum of remaining phyla from unfiltered provider abundances. Recorded finite precision may yield totals slightly different from 100%; values are not reclosed.","Descriptive only."),
        "Saved abundance-ranked taxa; overall and four-group means; no inference",width=180,height=100)
    })
  })
  plan("Fig3C",function(){
    z<-pub_species_upset(b,config)
    safe_write_csv(z$intersections,file.path(export,"plot_data/Fig3C_intersection_counts.csv"))
    safe_write_csv(z$set_sizes,file.path(export,"plot_data/Fig3C_set_sizes.csv"))
    save("Fig3C_BMI25_species_UpSet","Figure 3","C","Species overlap across five comparisons",z$plot,z$data,
      paste(basecaption,"Sets use MaAsLin2 full-model qval <0.05; intersections are exact memberships, not target-term BH."),
      "Exact set memberships; counts 69/32/10/4/5; intersection counts partition the union",q="qval",model="BMI25 species MaAsLin2",width=180,height=125)
  })
  plan("Fig3D",function(){
    d<-b$species[b$species$comparison==b$spec$comparison[1]&b$species$qval<.05,]
    d<-head(d[order(-abs(d$coef),d$qval,d$feature),],20)
    map<-read("04b_abundant_taxa/taxonomy_mapping.csv");map<-unique(map[!is.na(map$V3_short_name),c("V3_short_name","phylum")])
    stopifnot(!anyDuplicated(map$V3_short_name));d$phylum<-map$phylum[match(d$feature,map$V3_short_name)]
    stopifnot(!anyNA(d$phylum));d$label<-factor(pub_name(d$feature),levels=rev(pub_name(d$feature)))
    d$direction<-ifelse(d$coef>0,"Positive","Negative")
    p<-ggplot2::ggplot(d,ggplot2::aes(coef,label))+ggplot2::geom_segment(ggplot2::aes(x=0,xend=coef,yend=label),colour="grey75",linewidth=.4)+
      ggplot2::geom_point(ggplot2::aes(colour=direction,shape=phylum),size=2)+ggplot2::scale_colour_manual(values=pub_palette(cols,d$direction,"Fig3D direction"))+
      ggplot2::labs(x="T2D versus control coefficient",y=NULL,shape="Phylum",colour="Direction")+theme+ggplot2::theme(axis.text.y=ggplot2::element_text(face="italic"))
    save("Fig3D_BMI25_top20_absolute_effect","Figure 3","D","Largest species effects among significant features",p,d,
      "Twenty features selected among 69 full-model qval <0.05 associations by absolute coefficient; ties by qval and feature ID. This is a different selection from the clinical heatmap.",
      "Formal qval<0.05; descending abs(coef), ascending qval, feature ID; first20",q="qval",model="BMI25 overall species MaAsLin2",width=150,height=140)
  })
  plan("Fig3E",function(){
    ids<-sort(unique(b$species$feature[b$species$qval<.05]));grid<-expand.grid(feature=ids,comparison=b$spec$comparison)
    j<-match(paste(grid$feature,grid$comparison),paste(b$species$feature,b$species$comparison))
    grid$coef<-b$species$coef[j];grid$formal_q<-b$species$qval[j];grid$tested<-!is.na(j)
    grid$label<-factor(pub_name(grid$feature),levels=rev(pub_name(ids)));grid$comparison<-factor(grid$comparison,levels=b$spec$comparison)
    p<-pub_heatmap(grid,config,italic=TRUE)
    save("Fig3E_BMI25_species_effect_heatmap","Figure 3","E","Species associations across five contrasts",p,grid,
      paste(basecaption,"Union of full-model qval <0.05 features; asterisks indicate formal q values. Untested effects remain NA (grey), not zero."),
      "Union of all five formal significant sets; alphabetical feature IDs",q="qval",model="BMI25 species MaAsLin2",width=180,height=max(120,4.2*length(ids)))
  })
  for(i in 1:5)local({i<-i;panel<-c("A","B","C","E","D")[i]
    plan(paste0("Fig4",panel),function(){
      d<-b$species[b$species$comparison==b$spec$comparison[i],]
      d<-pub_select_volcano_labels(d,config$publication$volcano_labels$species_maximum,config$publication$volcano_labels$species_per_direction)
      save(paste0("Fig4",panel,"_BMI25_species_volcano"),"Figure 4",panel,b$spec$comparison[i],pub_volcano(d,config,d[d$label_selected,,drop=FALSE]),d,
        paste(basecaption,"MaAsLin2 full-model q value; all tested target coefficients are shown. Dashed line: qval=0.05."),
        "All tested coefficients; labels only q<0.05, at most five per direction and ten total; ascending q, descending abs(coef), feature ID",q="qval",model=paste(unique(d$model_id),collapse=";"),width=110,height=100)
    })
  })
  plan("Fig4F",function(){
    d<-b$clinical;sel<-b$clinical_selected[order(b$clinical_selected$feature_rank),]
    d$label<-factor(d$display_name,levels=rev(sel$display_name));d$trait<-factor(c(BMI="BMI",age="Age",sex="Sex")[d$clinical_trait],levels=c("BMI","Age","Sex"))
    d$stars<-pub_stars(d$q_bh_60);limit<-max(.25,abs(d$spearman_rho))
    p<-ggplot2::ggplot(d,ggplot2::aes(trait,label,fill=spearman_rho))+ggplot2::geom_tile(colour="white",linewidth=.2)+
      ggplot2::geom_text(ggplot2::aes(label=stars),size=2.8)+ggplot2::scale_fill_gradient2(low=cols["Negative"],mid="white",high=cols["Positive"],limits=c(-limit,limit),name="Spearman rho")+
      ggplot2::labs(x=NULL,y=NULL)+theme+ggplot2::theme(axis.text.y=ggplot2::element_text(face="italic"),axis.line=ggplot2::element_blank())
    save("Fig4F_Spearman_clinical_heatmap","Figure 4","F","Exploratory clinical correlations",p,d,
      "Twenty features ranked by overall MaAsLin2 full-model qval, pval and feature ID. Spearman correlations use relative abundance without additional transformation; Female=0 and Male=1. No covariate adjustment. Stars: BH-adjusted q value across 60 displayed tests, *<0.05, **<0.01, ***<0.001. Marginal associations do not establish independence or causality.",
      "Overall full species result sorted qval, pval, feature ID; first20; BMI, Age, Sex",q="q_bh_60",model="Exploratory Spearman 60-test family",width=120,height=145)
  })
  plan("Fig5A",function(){
    d<-pub_select_volcano_labels(b$pooled,config$publication$volcano_labels$pathway_maximum,
      config$publication$volcano_labels$pathway_per_direction,pathway=TRUE,root=run$root)
    save("Fig5A_MaAsLin2_overall_pathway_volcano","Figure 5","A","Pooled multivariable pathway associations",pub_volcano(d,config,d[d$label_selected,,drop=FALSE],parse_labels=FALSE),d,
      "All 5,945 pathway features; 1,641 associations with q<0.05. Fixed effects: disease, BMI, age, sex and study_name. TSS, LOG, standardize TRUE. Significance was based on the MaAsLin2 full-model q value.",
      "All 5945 target disease coefficients; labels only q<0.05, at most six per direction and twelve total; ascending q, descending abs(coef), feature ID; pathway ID plus abbreviated description",q="qval",model="Pooled pathway MaAsLin2",width=150,height=110)
  })
  plan("Fig5B",function(){
    d<-b$candidates[b$candidates$selected_for_genus_panel,]
    stopifnot(nrow(d)>0);count<-stats::aggregate(feature~genus,d,length);names(count)[2]<-"n_features"
    p<-ggplot2::ggplot(count,ggplot2::aes(n_features,reorder(genus,n_features)))+
      ggplot2::geom_col(fill=cols["Positive"],width=.65)+ggplot2::labs(x="Number of displayed pathway features",y=NULL)+theme+ggplot2::theme(axis.text.y=ggplot2::element_text(face="italic"))
    save("Fig5B_MaAsLin2_genus_pathway_features","Figure 5","B","Representative genus-stratified features",p,d,
      "Genus-level summary of HUMAnN taxonomic assignments among displayed taxon-stratified pathway features, including species-stratified assignments. This is not a separately fitted genus model. Candidate criterion abs(coef)>1 and full-model qval<0.01 is for visualization; formal significance remains qval<0.05.",
      paste("Candidates: abs(coef)>1 & qval<0.01; qval, descending abs(coef), feature ID; genus limit",config$publication$pathway_candidate_limit$genus),
      q="qval",model="Pooled pathway MaAsLin2",width=130,height=max(85,nrow(count)*5))
  })
  plan("Fig5C",function(){
    d<-b$candidates[b$candidates$selected_for_species_panel,];stopifnot(nrow(d)>0)
    d$full_pathway_name<-d$display_name
    d$pathway_label<-pub_wrap(paste0(d$base_pathway_ID,": ",pub_abbreviate(d$pathway_description,25)),35)
    assert_true(!anyDuplicated(d[,c("pathway_label","taxonomic_contributor")]),"Ambiguous species/pathway plot labels")
    p<-ggplot2::ggplot(d,ggplot2::aes(taxonomic_contributor,pathway_label,colour=coef,size=abs(coef)))+
      ggplot2::geom_point(alpha=.85)+ggplot2::scale_colour_gradient2(low=cols["Negative"],mid="white",high=cols["Positive"],midpoint=0)+
      ggplot2::scale_size_continuous(range=c(1.5,4))+ggplot2::labs(x=NULL,y=NULL,colour="Coefficient",size="|Coefficient|")+theme+
      ggplot2::theme(axis.text.x=ggplot2::element_text(angle=50,hjust=1,face="italic"),
        axis.text.y=ggplot2::element_text(size=7,lineheight=1.05),plot.margin=ggplot2::margin(7,10,7,18))
    save("Fig5C_MaAsLin2_species_pathway_features","Figure 5","C","Representative species-stratified features",p,d,
      "Taxon-stratified pathway contributions reflect HUMAnN assignment and inferred functional potential, not measured metabolic activity. Candidate criterion abs(coef)>1 and full-model qval<0.01; full candidate and display lists accompany Table S7.",
      paste("Candidates ranked qval, descending abs(coef), feature ID; species limit",config$publication$pathway_candidate_limit$species),q="qval",model="Pooled pathway MaAsLin2",width=180,height=min(180,max(110,length(unique(d$pathway_label))*6)))
  })
  plan("Fig6A",function(){
    a<-b$reml$pathway[b$reml$comparison==b$spec$comparison[1]&b$reml$padj<.05];z<-b$reml$pathway[b$reml$comparison==b$spec$comparison[2]&b$reml$padj<.05]
    d<-data.frame(set=c("Overall only","Shared","NT2D vs NC only"),n=c(length(setdiff(a,z)),length(intersect(a,z)),length(setdiff(z,a))))
    assert_true(identical(as.integer(d$n),c(9L,4L,129L)),"REML overlap counts differ")
    z<-pub_overlap_venn(d$n,c("Overall T2D vs Control","NT2D vs NC"),config)
    save("Fig6A_REML_pathway_overlap","Figure 6","A","Overlap of REML pathway associations",z$plot,d,
      "Overall 13, NT2D versus NC 133, shared 4 pathway features; three shared features assigned to Blautia wexlerae and one to Hungatella hathewayi. Formal padj is BH within each successful pooled contrast.",
      "Complete significant sets; non-area-proportional two-set Venn diagram",q="padj",model="DESeq2 + REML",width=150,height=105)
  })
  plan("Fig6B",function(){
    d<-pub_reml_heatmap_data(b,run$root)
    save("Fig6B_REML_pathway_effect_heatmap","Figure 6","B","Cross-cohort and BMI-stratified pathway patterns",pub_heatmap(d,config),d,
      "Five cohort-DESeq2 plus REML comparisons. Color: combined log2 fold change. Stars: contrast-specific padj. Untested values remain NA (grey), not zero. Patterns across strata do not establish interaction.",
      "All 13 overall significant features; 12 NT2D-only and 10 T2D-internal features by padj, raw P, feature ID; deduplicate, at most 35 rows",q="padj",model="DESeq2 + REML",width=180,height=max(140,length(unique(d$pathway))*4.5))
  })
  for(panel in c("C1","C2"))local({panel<-panel
    plan(paste0("Fig6",panel),function(){
      z<-final_pathway_forest(b,panel)
      save(paste0("Fig6",panel,"_REML_representative_forest"),"Figure 6",panel,z$comparison,
        ggplot2::ggplot(z$values,ggplot2::aes(coef,study))+ggplot2::geom_point(),z$values,
        final_forest_rule(),final_forest_rule(),q="padj",model="DESeq2 + REML",width=180,height=80)
    })
  })
  plan("Fig6D",function(){
    selected<-b$spec$comparison[c(1L,2L,5L)]
    d<-b$reml_annotation[b$reml_annotation$padj<.05 & b$reml_annotation$comparison%in%selected,]
    d$direction<-ifelse(d$combined_log2FC>0,"Positive","Negative")
    count<-as.data.frame(table(comparison=factor(d$comparison,levels=selected),functional_category=d$functional_category,direction=d$direction))
    count<-count[count$Freq>0,];count$comparison<-factor(count$comparison,levels=selected)
    direction_counts<-table(factor(d$comparison,levels=selected),factor(d$direction,levels=c("Positive","Negative")))
    assert_true(identical(as.integer(direction_counts),c(9L,114L,0L,4L,19L,20L)),"REML directional counts differ")
    p<-ggplot2::ggplot(count,ggplot2::aes(Freq,functional_category,fill=direction))+ggplot2::geom_col(width=.65)+
      ggplot2::facet_wrap(~comparison,nrow=1)+
      ggplot2::scale_fill_manual(values=pub_palette(cols,count$direction,"REML feature directions"))+
      ggplot2::labs(x="Significant pathway features",y=NULL,fill="Direction")+theme+
      ggplot2::theme(axis.text.y=ggplot2::element_text(size=7))
    safe_write_csv(d,file.path(export,"plot_data/Fig6D_significant_feature_classification.csv"))
    save("Fig6D_REML_functional_categories","Figure 6","D","Functional categories and directions of pathway features",p,count,
      "Counts describe taxon-stratified pathway features and inferred functional potential. Positive/negative counts: overall 9/4, NT2D versus NC 114/19, OT2D versus NT2D 0/20. Fixed functional categories are independent of direction; unmatched features are Other/Unclassified.",
      "All padj<0.05 features in three nonempty REML comparisons; category crossed with coefficient direction",q="padj",model="DESeq2 + REML",width=180,height=115)
  })
  plan("Fig7A",function(){
    d<-pub_interaction_summary(b)
    p<-ggplot2::ggplot(d,ggplot2::aes(domain_display,BMI_display))+ggplot2::geom_tile(fill="#E8EEF4",colour="white",linewidth=.8)+
      ggplot2::geom_text(ggplot2::aes(label=label),size=3,family=config$output$font_family)+ggplot2::labs(x=NULL,y=NULL)+theme+ggplot2::theme(axis.line=ggplot2::element_blank(),axis.ticks=ggplot2::element_blank())
    save("Fig7A_formal_interaction_families","Figure 7","A","Formal T2D by BMI interaction tests",p,d,
      "Number of FDR-significant target interaction terms / tested features. BH is separate for BMI25, continuous BMI and BMI23, and separately for species and pathway features. No FDR-significant feature interaction was detected.",
      "All six saved interaction families",q="formal_q",model="Formal T2D x BMI models",width=100,height=70)
  })
  plan("Fig7B",function(){
    d<-b$additive[b$additive$q_within_exposure_BH_BMI_per5<.05,];d<-d[order(d$coef_BMI_per5,d$feature),]
    d$label<-factor(pub_name(d$feature),levels=pub_name(d$feature));d$class<-ifelse(d$q_within_exposure_BH_T2D<.05,"Both","BMI-only")
    p<-ggplot2::ggplot(d,ggplot2::aes(coef_BMI_per5,label,colour=class))+
      ggplot2::geom_vline(xintercept=0,linetype=2,linewidth=.3)+ggplot2::geom_segment(ggplot2::aes(x=ci_lower_BMI_per5,xend=ci_upper_BMI_per5,yend=label),linewidth=.4)+
      ggplot2::geom_point(size=1.8)+ggplot2::scale_colour_manual(values=pub_palette(c(Both="#8961A3",`BMI-only`=cols[["Positive"]]),d$class,"Fig7B association class"))+
      ggplot2::labs(x="BMI coefficient per 5 kg/m² (95% CI)",y=NULL,colour="Association")+theme+ggplot2::theme(axis.text.y=ggplot2::element_text(face="italic"))
    save("Fig7B_additive_BMI_forest","Figure 7","B","BMI-associated species in the additive model",p,d,
      "Twenty-seven BMI-associated features in the additive model; exposure-specific BH q<0.05. Points and saved 95% confidence limits are per 5 kg/m². Both: also associated with T2D in its separate BH family.",
      "All 27 BMI term-specific BH-significant species; coefficient then feature ID",q="q_within_exposure_BH_BMI_per5",model="Additive T2D + BMI model",width=150,height=175)
  })
  plan("Fig7C",function(){
    d<-b$additive;t<-d$q_within_exposure_BH_T2D<.05;m<-d$q_within_exposure_BH_BMI_per5<.05
    d$class<-ifelse(t & m,"Both",ifelse(t,"T2D-only",ifelse(m,"BMI-only","Neither")))
    palette<-c(`T2D-only`=cols[["Negative"]],`BMI-only`=cols[["Positive"]],Both="#8961A3",Neither=cols[["NS"]])
    p<-ggplot2::ggplot(d,ggplot2::aes(coef_T2D,coef_BMI_per5,colour=class))+
      ggplot2::geom_hline(yintercept=0,colour="grey70",linewidth=.3)+ggplot2::geom_vline(xintercept=0,colour="grey70",linewidth=.3)+
      ggplot2::geom_point(size=1.5)+ggplot2::scale_colour_manual(values=pub_palette(palette,d$class,"Fig7C association class"),drop=FALSE)+ggplot2::labs(x="T2D coefficient",y="BMI coefficient per 5 kg/m²",colour="Association")+theme
    save("Fig7C_additive_effect_scatter","Figure 7","C","T2D and BMI additive associations",p,d,
      "All 200 species; T2D 62, BMI 27, both 13 and BMI-only 14. Categories use two exposure-specific BH families, not a joint 400-term family.",
      "All 200 species without label selection",q="Separate exposure-specific BH",model="Additive T2D + BMI model",width=120,height=100)
  })
  plan("Fig7D",function(){
    z<-pub_graphical_summary(b,config)
    save("Fig7D_association_summary","Figure 7","D","Interpretation of BMI-stratified associations",z$plot,z$data,
      "Differences in stratum-specific significance do not establish effect modification. Additive associations and FDR-negative interaction tests address different questions.",
      "Saved stratified counts, six formal interaction families and additive-set memberships",model="Interaction and additive summaries",width=180,height=125)
  })

  for(total in c(TRUE,FALSE))local({total<-total;panel<-if(total)"F" else "G"
    plan(paste0("SFig1",panel),function(){
      z<-pub_pcoa(read("04_primary_diversity/BMI25_beta_PCoA_coordinates.csv"),read("04_primary_diversity/BMI25_beta_sample_manifest.csv"),config,total)
      save(paste0("SFig1",panel,"_BMI25_ComBat_descriptive"),"Supplementary Figure 1",panel,"Descriptive community ordination",z$plot,z$data,
        "ComBat-corrected ordination is presented for descriptive visualization only. Dashed ellipses are 95% t data ellipses. Formal inference uses uncorrected Aitchison distances.",
        "Saved ComBat ordination coordinates; no differential-abundance use",model="Descriptive ComBat PCoA",width=110,height=85,supp=TRUE)
    })
  })
  for(rank in c("phylum","species"))local({rank<-rank;panel<-if(rank=="phylum")"A" else "B"
    plan(paste0("SFig2",panel),function(){
      d<-abundant[abundant$taxonomic_level==rank,,drop=FALSE]
      # A single implementation: save() routes this data into plot_abundant_taxa().
      # No inline legacy plot, model, emmeans call or multiple-testing correction.
      save(paste0("SFig2",panel,"_BMI25_top_",rank,"_boxplots"),"Supplementary Figure 2",panel,paste("Abundant",rank,"distributions"),
        NULL,d,paste(basecaption,"Secondary/exploratory abundant-taxa LMM; saved p_adj_BH <0.05 annotations only. Formal species inference remains MaAsLin2."),
        "Saved top4 phyla/top8 species; level-wide BH across 24/48 saved comparisons; no refit or adjustment",q="p_adj_BH",model="Abundant-taxa exploratory LMM",supp=TRUE)
    })
  })
  plan("SFig3",function(){
    z<-pub_sankey(b$candidates,config)
    save("SFig3_taxon_pathway_assignment_Sankey","Supplementary Figure 3","","Taxon-stratified pathway assignments",z$plot,z$data,
      "Links represent HUMAnN taxonomic assignment of selected pooled MaAsLin2 pathway features. Width counts displayed assignments; it does not represent metabolic flux or an experimentally verified mechanism. Candidate and display selection match Figure 5C.",
      "Same selected species-stratified candidate features as Figure5C",q="qval",model="Pooled pathway MaAsLin2 annotation",width=180,height=max(140,6*nrow(z$data)),supp=TRUE)
  })
  # Supplementary Figure4: three descriptive, all-matched-feature scatters.
  # Reuse the saved-model pairing and descriptive summaries used by S12/S13.
  # No model fit, q adjustment or correlation significance test is added here.
  sens<-publication_sensitivity_pairs(b,run)
  sensitivity_specs<-list(
    A=list(method="BMI23",domain="species",id="SFig4A_sensitivity_BMI23_species",
      title="BMI23 species coefficient concordance",xlab="BMI25 main-model coefficient",ylab="BMI23 sensitivity coefficient"),
    B=list(method="CLR",domain="species",id="SFig4B_sensitivity_CLR_species",
      title="CLR species coefficient concordance",xlab="BMI25 main-model coefficient",ylab="CLR species coefficient"),
    C=list(method="CLR",domain="pathways",id="SFig4C_sensitivity_CLR_pathways",
      title="CLR pathway coefficient concordance",xlab="BMI25 REML pooled log2 fold change",ylab="CLR pathway coefficient"))
  for(panel in names(sensitivity_specs))local({
    panel<-panel;spec<-sensitivity_specs[[panel]]
    plan(paste0("SFig4",panel),function(){
      d<-sens$pairs[sens$pairs$method==spec$method & sens$pairs$domain==spec$domain,,drop=FALSE]
      assert_true(nrow(d)>0L,paste("No sensitivity data:",spec$id))
      assert_true(!anyDuplicated(d[c("comparison","feature")]),"Duplicate matched sensitivity feature")
      assert_true(!anyNA(d$direction_agrees)&&all(is.finite(d$coef_main))&&all(is.finite(d$coef_sensitivity)),
        "Missing/nonfinite sensitivity coefficients; do not drop or relabel features")
      original<-d
      d$direction_class<-factor(ifelse(d$direction_agrees,"Concordant","Discordant"),levels=c("Concordant","Discordant"))
      d$shared_fdr<-!is.na(d$q_main)&!is.na(d$q_sensitivity)&d$q_main<.05&d$q_sensitivity<.05
      d$shared_status<-factor(ifelse(d$shared_fdr,"Shared FDR-significant","Not shared FDR-significant"),
        levels=c("Not shared FDR-significant","Shared FDR-significant"))
      # Keep source rows in exactly their previous order in plot-data. Draw pale
      # points first and dark points second using separate layers, not reselection.
      annotation<-sens$summary[sens$summary$method==spec$method & sens$summary$domain==spec$domain,,drop=FALSE]
      assert_true(!anyDuplicated(annotation$comparison)&&setequal(annotation$comparison,unique(d$comparison)),"Facet summary coverage mismatch")
      for(comparison in annotation$comparison){
        z<-d[d$comparison==comparison,,drop=FALSE];j<-match(comparison,annotation$comparison)
        assert_true(nrow(z)==annotation$n_matched[j]&&sum(z$shared_fdr)==annotation$shared_significant[j],"Saved descriptive counts differ")
      }
      r_display<-ifelse(is.finite(annotation$descriptive_coefficient_r),sprintf("%.2f",annotation$descriptive_coefficient_r),"NA")
      concordance_display<-ifelse(is.finite(annotation$direction_agreement),sprintf("%.1f%%",100*annotation$direction_agreement),"NA")
      annotation$label<-sprintf("n = %d\nConcordant = %s\nr = %s\nShared FDR n = %d",
        annotation$n_matched,concordance_display,r_display,annotation$shared_significant)
      direction_colours<-c(Concordant="#2166AC",Discordant="#F7AA58")
      shared_alpha<-c("Not shared FDR-significant"=.18,"Shared FDR-significant"=.95)
      plot<-final_sensitivity_scatter_plot(d,annotation,spec,config)
      assert_true(identical(d[names(original)],original),"Display preparation changed original paired values/order")
      safe_write_csv(annotation,file.path(export,"plot_data",paste0(spec$id,"_facet_summary.csv")))
      safe_write_csv(data.frame(direction_class=names(direction_colours),colour=unname(direction_colours)),
        file.path(export,"plot_data",paste0(spec$id,"_direction_palette.csv")))
      safe_write_csv(data.frame(shared_status=names(shared_alpha),alpha=unname(shared_alpha)),
        file.path(export,"plot_data",paste0(spec$id,"_shared_alpha.csv")))
      caption<-paste("Each point represents a feature matched by exact feature ID between the BMI25 main analysis and the sensitivity analysis.",
        "Blue points indicate concordant coefficient directions and orange points indicate discordant directions.",
        "Dark points were FDR significant in both analyses; pale points were not jointly FDR significant.",
        "Facet annotations report matched feature count, direction concordance, descriptive coefficient correlation r, and shared FDR-significant feature count.",
        "These summaries are descriptive and are not tests of equivalence or effect modification. No correlation P value is calculated.")
      if(spec$method=="CLR")caption<-paste(caption,"The main and CLR coefficients are on different scales and should not be interpreted as numerically equivalent.")
      save(spec$id,"Supplementary Figure 4",panel,spec$title,plot,d,caption,
        "All features matched by exact feature ID; original row order preserved; no top-N selection, no q recalculation and no new statistical tests",
        q="Saved main-model and sensitivity-model FDR columns",model=paste(spec$method,spec$domain,"sensitivity comparison"),width=config$publication$sensitivity_scatter$width_mm,height=config$publication$sensitivity_scatter$height_mm,supp=TRUE)
    })
  })
  plan("SFig1C_combined",function(){
    z<-final_cohort_alpha(b,run,refinement)
    save("SFig1C_BMI25_cohort_alpha_combined","Supplementary Figure 1","C","Alpha Diversity by Cohort",z$plot,z$data,
      "Upper: Shannon diversity; lower: species richness. The same 1186 participants are displayed for each metric. Metric-specific pooled LMM P values and ICC are shown; no new tests are computed. A common group legend and one horizontal cohort axis are displayed.",
      "All current-run alpha observations; two metrics stacked with common legend",model="Alpha diversity LMM",width=z$width,height=z$height,supp=TRUE)
  })
  checks<-pub_bind(preflight)
  safe_write_csv(checks,file.path(export,"validation/panel_preparation.csv"))
  required<-pub_read_csv_snapshot(run,file.path(run$root,relative_project_path(config$publication$panel_registry)))
  pub_key_contract(required,"figure_id","Panel registry")
  pub_validate_panel_counts(required,config)
  missing<-setdiff(required$figure_id,names(pending));extra<-setdiff(names(pending),required$figure_id)
  safe_write_csv(data.frame(figure_id=required$figure_id,prepared=required$figure_id%in%names(pending)),
    file.path(export,"validation/panel_coverage.csv"))
  assert_true(!length(extra),paste("Unregistered panel IDs:",paste(extra,collapse=", ")))
  if(isTRUE(config$publication$require_all_panels))assert_true(all(checks$preparation_complete)&&!length(missing),
    paste("Panel preflight incomplete; no image devices opened. See validation/panel_preparation.csv and panel_coverage.csv. Missing:",paste(missing,collapse=", ")))
  for(id in names(pending))local({id<-id;attempt(id,"figure",function(){
    index[[id]]<<-do.call(pub_save_panel,pending[[id]])
  })})
  safe_write_csv(pub_bind(index),file.path(export,"figure_index.csv"))
  pub_finalize_figure_integrity(export)
  coverage<-data.frame(figure_id=required$figure_id,exported=required$figure_id%in%names(index))
  safe_write_csv(coverage,file.path(export,"validation/panel_export_coverage.csv"))
  if(isTRUE(config$publication$require_all_panels))assert_true(all(coverage$exported),
    "Some planned panels were not exported; inspect panel_export_coverage.csv and export_gaps.csv")
  pub_validate_panel_counts(pub_bind(index),config)
  pub_write_combined_legends(pub_bind(index),export)
  invisible(index)
}
