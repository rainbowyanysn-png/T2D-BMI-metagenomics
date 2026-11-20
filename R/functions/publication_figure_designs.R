# Presentation-only selections and geometry; statistics always come from saved results.
pub_select_volcano_labels <- function(x, maximum = 10L, per_direction = 5L, pathway = FALSE, root = NULL) {
  pub_key_contract(x, "feature", "Volcano label inputs")
  ix <- which(!is.na(x$qval) & x$qval < .05 & is.finite(x$coef))
  ix <- ix[order(x$qval[ix], -abs(x$coef[ix]), x$feature[ix])]
  selected <- c(head(ix[x$coef[ix] > 0], per_direction), head(ix[x$coef[ix] < 0], per_direction))
  selected <- head(selected[order(x$qval[selected], -abs(x$coef[selected]), x$feature[selected])], maximum)
  x$label_selected <- seq_len(nrow(x)) %in% selected
  x$label_rank <- match(seq_len(nrow(x)), selected)
  x$label_text <- if (pathway) {
    annotation <- annotate_pathway_features(x$pathway, root)
    paste0(annotation$base_pathway_ID, ": ", pub_abbreviate(annotation$pathway_description, 27L))
  } else species_plot_label(x$feature)
  assert_true(sum(x$label_selected) <= maximum && all(x$qval[x$label_selected] < .05), "Invalid volcano labels")
  x
}

pub_abbreviate <- function(x, width = 32L) {
  x <- gsub("[._]+", " ", x)
  ifelse(nchar(x) > width, paste0(substr(x, 1L, width - 3L), "..."), x)
}

pub_workflow_design <- function(b, config) {
  d <- data.frame(module = c("Data sources", "Study population", "Data processing", "Primary analyses",
    "Additional inference", "Sensitivity analyses"), x = rep(1:3, 2), y = rep(c(3, 1), each = 3),
    detail = c("Five curatedMetagenomicData\ncohorts",
      "1,186 participants\nNC 307   OC 119\nNT2D 138   OT2D 622",
      "Species and MetaCyc pathway features\nBMI25 groups\nPrespecified prevalence filters",
      "Alpha and beta diversity\nFive MaAsLin2 species contrasts\nPooled MaAsLin2 pathway features\nCohort DESeq2 + REML",
      "Formal T2D-by-BMI interactions\nT2D + BMI additive model\nExploratory Spearman correlations",
      "Cohort-adapted BMI23\nCZM + CLR\nMedication-adjustment sensitivity"))
  edges <- data.frame(x = c(1.43, 2.43, 3, 1, 2, 3), y = c(3, 3, 2.49, 2, 2, 2),
    xend = c(1.57, 2.57, 3, 1, 2, 3), yend = c(3, 3, 2, 1.57, 1.57, 1.57))
  p <- ggplot2::ggplot(d) +
    ggplot2::geom_segment(data = edges, ggplot2::aes(x=x,y=y,xend=xend,yend=yend),
      arrow=grid::arrow(length=grid::unit(1.3,"mm")), colour="#376795", linewidth=.45) +
    ggplot2::annotate("segment", x=1,xend=3,y=2,yend=2,colour="#376795",linewidth=.45) +
    ggplot2::geom_rect(ggplot2::aes(xmin=x-.43,xmax=x+.43,ymin=y-.5,ymax=y+.5),
      fill="#EAF1F6",colour="#376795",linewidth=.45) +
    ggplot2::geom_text(ggplot2::aes(x,y=y+.34,label=module),fontface="bold",size=2.9,family=config$output$font_family) +
    ggplot2::geom_text(ggplot2::aes(x,y=y-.04,label=detail),size=2.5,lineheight=1.2,family=config$output$font_family) +
    ggplot2::coord_cartesian(xlim=c(.5,3.5),ylim=c(.4,3.6),clip="off") +
    ggplot2::theme_void(base_family=config$output$font_family)
  list(plot=p,data=d)
}

pub_species_upset <- function(b, config) {
  features <- sort(unique(b$species$feature[b$species$qval < .05]))
  membership <- sapply(b$spec$comparison, function(cmp)
    features %in% b$species$feature[b$species$comparison == cmp & b$species$qval < .05])
  assert_true(identical(as.integer(colSums(membership)),c(69L,32L,10L,4L,5L)), "UpSet set sizes differ")
  keys <- apply(membership,1,paste,collapse="|")
  counts <- table(keys); ord <- order(-as.integer(counts),names(counts)); counts <- counts[ord]
  assert_true(sum(counts)==length(features),"UpSet intersections do not partition the union")
  bar <- data.frame(intersection=seq_along(counts),n=as.integer(counts),pattern=names(counts))
  dots <- pub_bind(lapply(seq_len(nrow(bar)),function(i)data.frame(intersection=i,
    comparison=factor(b$spec$comparison,levels=rev(b$spec$comparison)),
    member=strsplit(bar$pattern[i],"|",fixed=TRUE)[[1]]=="TRUE")))
  links <- pub_bind(lapply(seq_len(nrow(bar)),function(i){z<-dots[dots$intersection==i & dots$member,,drop=FALSE];
    data.frame(intersection=i,lower=min(as.integer(z$comparison)),upper=max(as.integer(z$comparison)))}))
  sizes <- data.frame(comparison=factor(b$spec$comparison,levels=rev(b$spec$comparison)),n=as.integer(colSums(membership)))
  theme <- publication_panel_theme(config)
  top <- ggplot2::ggplot(bar,ggplot2::aes(intersection,n))+ggplot2::geom_col(fill="#2166AC",width=.7)+
    ggplot2::geom_text(ggplot2::aes(label=n),vjust=-.35,size=2.5)+
    ggplot2::scale_x_continuous(limits=c(.4,nrow(bar)+.6),expand=c(0,0))+
    ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(0,.15)))+
    ggplot2::labs(x=NULL,y="Intersection size")+theme+
    ggplot2::theme(axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank())
  bottom <- ggplot2::ggplot(dots,ggplot2::aes(intersection,as.integer(comparison)))+
    ggplot2::geom_segment(data=links,ggplot2::aes(x=intersection,xend=intersection,y=lower,yend=upper),inherit.aes=FALSE,linewidth=.4)+
    ggplot2::geom_point(ggplot2::aes(colour=member),size=1.9)+
    ggplot2::scale_colour_manual(values=c(`TRUE`="#222222",`FALSE`="#DDDDDD"))+
    ggplot2::scale_x_continuous(limits=c(.4,nrow(bar)+.6),expand=c(0,0))+
    ggplot2::scale_y_continuous(limits=c(.5,5.5),breaks=1:5,expand=c(0,0))+
    ggplot2::labs(x="Distinct intersections",y=NULL)+theme+
    ggplot2::theme(legend.position="none",axis.text.y=ggplot2::element_blank(),axis.ticks=ggplot2::element_blank(),axis.line=ggplot2::element_blank())
  left <- ggplot2::ggplot(sizes,ggplot2::aes(n,comparison))+ggplot2::geom_col(fill="#528FAD",width=.6)+
    ggplot2::geom_text(ggplot2::aes(label=n),hjust=-.15,size=2.4)+
    ggplot2::scale_x_reverse(expand=ggplot2::expansion(mult=c(.25,0)))+
    ggplot2::scale_y_discrete(expand=ggplot2::expansion(add=.5))+
    ggplot2::labs(x="Set size",y=NULL)+theme+ggplot2::theme(axis.text.y=ggplot2::element_text(size=7))
  p <- patchwork::wrap_plots(patchwork::plot_spacer(),top,left,bottom,ncol=2,widths=c(1.8,3.8),heights=c(2,1.3))
  raw <- data.frame(feature=features,membership,intersection_pattern=keys,
    intersection_count=as.integer(counts[match(keys,names(counts))]),check.names=FALSE)
  list(plot=p,data=raw,intersections=bar,set_sizes=sizes)
}

pub_overlap_venn <- function(counts, labels, config, colours=c("#72BCD5","#F7AA58")) {
  assert_true(length(counts)==3L && all(counts>=0),"Invalid two-set overlap counts")
  t <- seq(0,2*pi,length.out=240)
  circles <- pub_bind(lapply(1:2,function(i)data.frame(x=c(-.55,.55)[i]+cos(t),y=sin(t),set=labels[i])))
  text <- data.frame(x=c(-1.02,0,1.02),y=0,label=as.character(counts))
  p <- ggplot2::ggplot()+ggplot2::geom_polygon(data=circles,ggplot2::aes(x,y,group=set,fill=set),alpha=.48,colour="white",linewidth=.4)+
    ggplot2::geom_text(data=text,ggplot2::aes(x,y,label=label),size=4,family=config$output$font_family)+
    ggplot2::annotate("text",x=c(-.75,.75),y=1.22,label=labels,size=2.7,family=config$output$font_family)+
    ggplot2::scale_fill_manual(values=pub_palette(setNames(colours,labels),circles$set,"Overlap sets"))+
    ggplot2::coord_equal(xlim=c(-1.65,1.65),ylim=c(-1.1,1.5),clip="off")+
    ggplot2::theme_void(base_family=config$output$font_family)+ggplot2::theme(legend.position="none")
  list(plot=p,data=data.frame(region=c("First only","Shared","Second only"),n=counts,first=labels[1],second=labels[2]))
}

pub_reml_heatmap_data <- function(b, root, full=FALSE) {
  r <- b$reml
  ranked <- function(cmp) {z<-r[r$comparison==cmp & r$padj<.05,,drop=FALSE];z[order(z$padj,z$p_value,z$pathway),,drop=FALSE]}
  overall <- ranked(b$spec$comparison[1]); nt <- ranked(b$spec$comparison[2]); internal <- ranked(b$spec$comparison[5])
  selected <- list(overall=overall$pathway,NT2D_only=head(nt$pathway[!nt$pathway%in%overall$pathway],12L),
    T2D_internal=head(internal$pathway,10L))
  if(full){z<-r[r$padj<.05,,drop=FALSE];z<-z[order(z$padj,z$p_value,z$pathway),];ids<-unique(z$pathway)} else {
    ids<-unique(unlist(selected,use.names=FALSE));assert_true(length(ids)<=35L && all(overall$pathway%in%ids),"Invalid representative REML heatmap selection")
  }
  d<-expand.grid(pathway=ids,comparison=b$spec$comparison,stringsAsFactors=FALSE)
  j<-match(paste(d$pathway,d$comparison),paste(r$pathway,r$comparison))
  d$coef<-r$combined_log2FC[j];d$formal_q<-r$padj[j];d$tested<-!is.na(j)
  ann<-annotate_pathway_features(ids,root)
  short<-paste(ann$base_pathway_ID,ifelse(is.na(ann$taxonomic_contributor),"unstratified",ann$taxonomic_contributor),sep=" | ")
  assert_true(!anyDuplicated(short),"Ambiguous REML heatmap row labels")
  d$label<-factor(short[match(d$pathway,ids)],levels=rev(short))
  d$comparison<-factor(d$comparison,levels=b$spec$comparison)
  d$selection_overall<-d$pathway%in%selected$overall
  d$selection_NT2D_only<-d$pathway%in%selected$NT2D_only
  d$selection_T2D_internal<-d$pathway%in%selected$T2D_internal
  d$full_pathway_name<-d$pathway
  assert_true(all(is.na(d$coef[!d$tested])) && all(is.na(d$formal_q[!d$tested])),"Untested heatmap values must remain missing")
  d
}

pub_interaction_summary <- function(b) {
  d<-pub_bind(lapply(split(b$interactions,list(b$interactions$BMI_definition,b$interactions$domain),drop=TRUE),function(z)
    data.frame(BMI_definition=z$BMI_definition[1],domain=z$domain[1],n_tested=nrow(z),n_significant=sum(z$formal_q<.05))))
  d$label<-paste0(d$n_significant,"/",format(d$n_tested,big.mark=",",trim=TRUE))
  d$BMI_display<-factor(c(BMI25="BMI25",continuous_BMI="Continuous BMI",BMI23="Cohort-adapted BMI23")[d$BMI_definition],
    levels=c("Cohort-adapted BMI23","Continuous BMI","BMI25"))
  d$domain_display<-factor(c(species="Species",pathways="Pathway features")[d$domain],levels=c("Species","Pathway features"))
  assert_true(nrow(d)==6L && all(d$n_significant==0L),"Interaction summary must retain six FDR-negative families")
  d
}

pub_graphical_summary <- function(b,config) {
  labels<-b$spec$comparison[c(2,3,5)]
  left<-pub_bind(lapply(labels,function(cmp)data.frame(comparison=cmp,domain=c("Species","Pathway features"),
    n=c(sum(b$species$qval[b$species$comparison==cmp]<.05),sum(b$reml$padj[b$reml$comparison==cmp]<.05)))))
  left$comparison<-factor(left$comparison,levels=rev(labels))
  theme<-publication_panel_theme(config)
  p1<-ggplot2::ggplot(left,ggplot2::aes(n,comparison,fill=domain))+
    ggplot2::geom_col(position=ggplot2::position_dodge(width=.7),width=.6)+
    ggplot2::geom_text(ggplot2::aes(label=n),position=ggplot2::position_dodge(width=.7),hjust=-.2,size=2.5)+
    ggplot2::scale_fill_manual(values=c(Species="#528FAD",`Pathway features`="#F7AA58"))+
    ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=c(0,.2)))+
    ggplot2::labs(x="Associated features",y=NULL,title="BMI-stratified\nassociation patterns",fill=NULL)+theme+
    ggplot2::theme(legend.position="bottom",plot.title=ggplot2::element_text(size=8))
  middle<-pub_interaction_summary(b)
  p2<-ggplot2::ggplot(middle,ggplot2::aes(domain_display,BMI_display))+
    ggplot2::geom_tile(fill="#EAF1F6",colour="white",linewidth=1)+
    ggplot2::geom_text(ggplot2::aes(label=label),size=2.7,family=config$output$font_family)+
    ggplot2::labs(x=NULL,y=NULL,title="No FDR-significant\nT2D-by-BMI interaction")+theme+
    ggplot2::theme(axis.line=ggplot2::element_blank(),axis.ticks=ggplot2::element_blank(),axis.text.x=ggplot2::element_text(angle=25,hjust=1),plot.title=ggplot2::element_text(size=8))
  a<-b$additive;t<-a$q_within_exposure_BH_T2D<.05;m<-a$q_within_exposure_BH_BMI_per5<.05
  counts<-c(sum(t & !m),sum(t & m),sum(!t & m))
  right<-pub_overlap_venn(counts,c("T2D (62)","BMI (27)"),config)
  p3<-right$plot+ggplot2::labs(title="Partly overlapping\nadditive associations")+
    ggplot2::theme(plot.title=ggplot2::element_text(size=8,hjust=.5,family=config$output$font_family))
  conclusion<-"BMI provides complementary phenotypic context, while statistical effect modification was not established."
  p<-patchwork::wrap_plots(p1,p2,p3,nrow=1,widths=c(1.2,1.25,1))+
    patchwork::plot_annotation(caption=conclusion,theme=ggplot2::theme(plot.caption=ggplot2::element_text(size=8,hjust=.5,family=config$output$font_family)))
  data<-pub_bind(list(transform(left,section="Stratified associations"),
    transform(middle,section="Formal interactions"),transform(right$data,section="Additive associations")))
  list(plot=p,data=data)
}
