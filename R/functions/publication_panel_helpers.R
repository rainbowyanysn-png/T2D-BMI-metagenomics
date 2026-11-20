pub_effect_colours <- function(config) pub_palette(config$publication$effect_colours,c("Positive","Negative","NS"),"effect directions")
pub_stars <- function(q) ifelse(is.na(q),"",ifelse(q<.001,"***",ifelse(q<.01,"**",ifelse(q<.05,"*",""))))
pub_name <- function(x) clinical_feature_display(x)
pub_wrap <- function(x,width=42) vapply(x,function(z)paste(strwrap(z,width=width),collapse="\n"),character(1))

pub_save_panel <- function(plot,id,figure,panel,title,data,caption,rule,export,config,
                           q="Not applicable",model="Descriptive",width=NULL,height=NULL,supp=FALSE) {
  title<-pub_caption_text(title);caption<-pub_caption_text(caption);rule<-pub_caption_text(rule)
  pub_language_check(c(title,caption,rule),id)
  assert_true(length(id)==1L&&grepl("^[A-Za-z0-9_]+$",id),"Panel ID must be a simple registered filename stem")
  width <- if(is.null(width))config$publication$half_width_mm else width
  height <- if(is.null(height))config$publication$default_height_mm else height
  tiff_dpi<-config$publication$tiff_dpi_by_figure[[id]]
  if(is.null(tiff_dpi))tiff_dpi<-300
  assert_true(length(tiff_dpi)==1L&&is.finite(tiff_dpi)&&tiff_dpi>=300,"Invalid TIFF DPI setting")
  pub_validate_panel(plot,data,id,width,height)
  if("label_selected"%in%names(data)){
    drawn<-unlist(lapply(plot$layers,function(layer)if(inherits(layer$geom,"GeomTextRepel"))layer$data$feature else character()),use.names=FALSE)
    assert_true(identical(sort(drawn),sort(data$feature[data$label_selected])),"Drawn volcano labels differ from saved selection")
  }
  plot <- plot + ggplot2::labs(tag=NULL)
  folder <- if(supp)"supplementary_panels" else "main_panels"
  paths <- file.path(export,"figures",folder,paste0(id,c(".pdf",".tiff")))
  png <- file.path(export,"figures/previews",paste0(id,".png"))
  stopifnot(!any(file.exists(c(paths,png))))
  staging<-file.path(export,"logs/panel_staging",id)
  assert_true(!dir.exists(staging),paste("Refusing to reuse panel staging:",id))
  dir.create(staging,recursive=TRUE)
  staged<-file.path(staging,paste0(id,c(".pdf",".tiff",".png")))
  ggplot2::ggsave(staged[1],plot,width=width,height=height,units="mm",bg="white",device=grDevices::cairo_pdf,limitsize=FALSE)
  ggplot2::ggsave(staged[2],plot,width=width,height=height,units="mm",dpi=tiff_dpi,bg="white",device=grDevices::tiff,type="cairo",compression="lzw",limitsize=FALSE)
  ggplot2::ggsave(staged[3],plot,width=width,height=height,units="mm",dpi=300,bg="white",device=grDevices::png,type="cairo",limitsize=FALSE)
  dat <- file.path(export,"plot_data",paste0(id,".csv"));cap <- file.path(export,"figures/captions",paste0(id,".txt"))
  rule_path<-file.path(export,"plot_data",paste0(id,"_selection_rule.txt"))
  metadata<-file.path(staging,c("plot_data.csv","caption.txt","selection_rule.txt"))
  safe_write_csv(data,metadata[1]);safe_write_lines(c(title,caption),metadata[2]);safe_write_lines(rule,metadata[3])
  from<-c(staged,metadata);to<-c(paths,png,dat,cap,rule_path)
  assert_true(!any(file.exists(to))&&all(file.exists(from))&&all(file.info(from)$size>0),paste("Incomplete staged panel:",id))
  for(i in seq_along(from))assert_true(file.rename(from[i],to[i]),paste("Cannot publish staged file:",to[i]))
  rel <- function(x)substring(x,nchar(export)+2L)
  record<-data.frame(figure_id=id,manuscript_figure=figure,panel=panel,title=title,PDF_path=rel(paths[1]),
    TIFF_path=rel(paths[2]),PNG_path=rel(png),plot_data_path=rel(dat),caption_path=rel(cap),
    width_mm=width,height_mm=height,dpi=tiff_dpi,legend_status="Included when applicable",
    formal_q_column=q,source_model_id=model,selection_rule=rule,selection_rule_path=rel(rule_path),
    tiff_compression="LZW",background="white")
  pub_validate_saved_panel(record,export)
  record
}

pub_volcano <- function(x,config,labels=NULL,parse_labels=TRUE) {
  x$state<-ifelse(x$qval<.05,ifelse(x$coef>0,"Positive","Negative"),"NS")
  x$log_q<- -log10(pmax(x$qval,.Machine$double.xmin))
  p<-ggplot2::ggplot(x,ggplot2::aes(coef,log_q,colour=state))+
    ggplot2::geom_point(size=.9,alpha=.75)+ggplot2::geom_hline(yintercept=-log10(.05),linetype=2,linewidth=.3)+
    ggplot2::geom_vline(xintercept=0,colour="grey60",linewidth=.25)+
    ggplot2::scale_colour_manual(values=pub_palette(pub_effect_colours(config),x$state,"volcano direction"),drop=FALSE)+
    ggplot2::labs(x="MaAsLin2 coefficient",y=expression(-log[10]("q value")),colour="Association")+publication_panel_theme(config)
  if(!is.null(labels) && nrow(labels)) {
    assert_true(all(labels$qval<.05),"Volcano labels must use formal full-model q significance")
    labels$log_q<- -log10(pmax(labels$qval,.Machine$double.xmin))
    p<-p+ggrepel::geom_text_repel(data=labels,ggplot2::aes(coef,log_q,label=label_text),
      inherit.aes=FALSE,parse=parse_labels,seed=config$publication$volcano_labels$seed,
      family=config$output$font_family,size=2.4,max.overlaps=Inf,min.segment.length=0,
      box.padding=.35,point.padding=.2,segment.size=.25)
  }
  p
}

pub_heatmap <- function(x,config,effect="coef",q="formal_q",label="label",xcol="comparison",italic=FALSE) {
  x$effect_plot<-x[[effect]];x$q_label<-pub_stars(x[[q]])
  limit<-max(abs(x$effect_plot),na.rm=TRUE)
  if(!is.finite(limit)||limit==0)limit<-1
  p<-ggplot2::ggplot(x,ggplot2::aes(.data[[xcol]],.data[[label]],fill=effect_plot))+
    ggplot2::geom_tile(colour="white",linewidth=.15)+ggplot2::geom_text(ggplot2::aes(label=q_label),size=2.4)+
    ggplot2::scale_fill_gradient2(low=pub_effect_colours(config)["Negative"],mid="white",
      high=pub_effect_colours(config)["Positive"],limits=c(-limit,limit),na.value="#E0E0E0",name="Coefficient")+
    ggplot2::labs(x=NULL,y=NULL)+publication_panel_theme(config)+
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1),axis.line=ggplot2::element_blank())
  if(italic)p<-p+ggplot2::theme(axis.text.y=ggplot2::element_text(face="italic"))
  p
}

pub_pcoa <- function(coords,md,config,total=FALSE) {
  stopifnot(!anyDuplicated(coords$sample_id),!anyDuplicated(md$sample_id))
  j<-match(coords$sample_id,md$sample_id);stopifnot(!anyNA(j))
  d<-coords;d$study_name<-md$study_name[j]
  d$phenotype<-if(total)factor(ifelse(md$disease[j]=="control","Control","T2D"),levels=c("Control","T2D")) else
    factor(unname(group_display[md$analysis_group[j]]),levels=c("NC","OC","NT2D","OT2D"))
  style<-read_presentation_config(config$project_root)
  palette<-unlist(if(total)style$groups$total else style$groups$fill)
  if(total)names(palette)[names(palette)=="Ctrl"]<-"Control"
  shapes<-setNames(as.integer(unlist(style$beta$cohort_shapes)),unlist(config$cohorts$names))
  p<-ggplot2::ggplot(d,ggplot2::aes(PC1,PC2,colour=phenotype))+
    ggplot2::geom_point(ggplot2::aes(shape=study_name),size=.9,alpha=.55)+
    ggplot2::stat_ellipse(ggplot2::aes(group=phenotype),type=style$beta$ellipse_type,level=style$beta$ellipse_level,
      linetype=2,linewidth=.4,show.legend=FALSE)+
    ggplot2::scale_colour_manual(values=pub_palette(palette,d$phenotype,"PCoA phenotype"),drop=FALSE)+ggplot2::scale_shape_manual(values=pub_named_mapping(shapes,d$study_name,"PCoA cohort shapes"))+
    ggplot2::labs(colour="Group",shape="Cohort",x="PCoA 1",y="PCoA 2")+publication_panel_theme(config)
  if(all(c("PC1_fraction","PC2_fraction")%in%names(d)))p<-p+ggplot2::labs(
    x=sprintf("PCoA 1 (%.1f%%)",100*d$PC1_fraction[1]),y=sprintf("PCoA 2 (%.1f%%)",100*d$PC2_fraction[1]))
  list(plot=p,data=d)
}

pub_diagnostic_plot <- function(x,config,kind="PERMANOVA") {
  pcol<-intersect(c("p_interaction","Pr(>F)","P_value","p_value","Pr(>Chisq)"),names(x))[1]
  value<-intersect(c("R2","F","F value","F_value","coefficient","P_value"),names(x))[1]
  stopifnot(!is.na(pcol),!is.na(value))
  labelcols<-intersect(c("metric","interaction_term","study_name","term","group_variable","Analysis","Variable"),names(x))
  stopifnot(length(labelcols)>0)
  x$test_label<-apply(x[,labelcols,drop=FALSE],1,paste,collapse=" / ")
  x$value_plot<-x[[value]];x$p_label<-paste0("P ",saved_p_label(x[[pcol]]))
  p<-ggplot2::ggplot(x,ggplot2::aes(value_plot,factor(test_label,levels=rev(unique(test_label)))))+
    ggplot2::geom_point(colour=pub_effect_colours(config)["Negative"],size=1.8)+
    ggplot2::geom_text(ggplot2::aes(label=p_label),hjust=-.1,size=2.5,family=config$output$font_family)+
    ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=c(.05,.65)))+
    ggplot2::labs(x=paste(kind,value),y=NULL)+publication_panel_theme(config)
  if(value=="coefficient" && all(c("ci_lower","ci_upper")%in%names(x)))p<-p+
    ggplot2::geom_segment(ggplot2::aes(x=ci_lower,xend=ci_upper,
      y=factor(test_label,levels=rev(unique(test_label))),yend=factor(test_label,levels=rev(unique(test_label)))),linewidth=.4)
  list(plot=p,data=x)
}

pub_sankey <- function(candidates,config) {
  d<-candidates[candidates$selected_for_species_panel,]
  stopifnot(nrow(d)>0L,!anyNA(d$genus),!anyNA(d$taxonomic_contributor))
  # Each selected feature contributes one assignment, not an abundance/flux weight.
  axes<-list(species=d$taxonomic_contributor,genus=d$genus,pathway=d$base_pathway_ID)
  nodes<-list();pos<-list()
  for(k in 1:3){
    count<-table(axes[[k]]);names_order<-sort(names(count));gap<-.018
    heights<-as.numeric(count[names_order])/sum(count)*(1-gap*(length(count)-1))
    lower<-c(0,head(cumsum(heights+gap),-1));upper<-lower+heights
    z<-data.frame(axis=k,node=names_order,ymin=lower,ymax=upper,xmin=k-.025,xmax=k+.025)
    nodes[[k]]<-z;pos[[k]]<-setNames(rep(0,length(names_order)),names_order)
  }
  ribbons<-list();n<-0L
  for(k in 1:2)for(i in seq_len(nrow(d))){
    a<-axes[[k]][i];z<-axes[[k+1]][i];left<-nodes[[k]][nodes[[k]]$node==a,];right<-nodes[[k+1]][nodes[[k+1]]$node==z,]
    lh<-(left$ymax-left$ymin)/sum(axes[[k]]==a);rh<-(right$ymax-right$ymin)/sum(axes[[k+1]]==z)
    # Left and right counters are separate for each pair of adjacent axes.
    if(i==1L){ lc<-setNames(rep(0,nrow(nodes[[k]])),nodes[[k]]$node);rc<-setNames(rep(0,nrow(nodes[[k+1]])),nodes[[k+1]]$node) }
    lo1<-left$ymin+lc[a]*lh;lo2<-right$ymin+rc[z]*rh;lc[a]<-lc[a]+1;rc[z]<-rc[z]+1
    t<-seq(0,1,length.out=40);smooth<-3*t^2-2*t^3
    n<-n+1L;ribbons[[n]]<-data.frame(x=c(k+.025+t*.95,rev(k+.025+t*.95)),
      y=c(lo1+(lo2-lo1)*smooth,rev(lo1+lh+(lo2+rh-lo1-lh)*smooth)),ribbon=n,genus=d$genus[i])
  }
  nd<-pub_bind(nodes);nd$face<-ifelse(nd$axis<3,"italic","plain")
  p<-ggplot2::ggplot()+ggplot2::geom_polygon(data=pub_bind(ribbons),ggplot2::aes(x,y,group=ribbon,fill=genus),alpha=.55,linewidth=0)+
    ggplot2::geom_rect(data=nd,ggplot2::aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="#8796A2",colour="white",linewidth=.2)+
    ggplot2::geom_text(data=nd,ggplot2::aes(axis,(ymin+ymax)/2,label=node,fontface=face),size=2.3,family=config$output$font_family)+
    ggplot2::scale_fill_manual(values=pub_palette(setNames(grDevices::hcl.colors(length(unique(d$genus)),"Dark 3"),sort(unique(d$genus))),d$genus,"Sankey genus assignment"))+
    ggplot2::scale_x_continuous(breaks=1:3,labels=c("Species","Genus","Pathway ID"),expand=ggplot2::expansion(mult=c(.12,.12)))+
    ggplot2::labs(x=NULL,y=NULL)+publication_panel_theme(config)+
    ggplot2::theme(axis.text.y=ggplot2::element_blank(),axis.ticks.y=ggplot2::element_blank(),axis.line=ggplot2::element_blank(),legend.position="none")
  list(plot=p,data=d)
}
