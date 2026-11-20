plot_composition <- function(d,panel,root,paths,style) {
 r2_columns(d,c("scope","group","feature","raw_abundance","display_name"))
 species<-panel$plot_family=="composition_species"
 r1_need(nrow(d)==if(species)54L else 30L,"Composition row count differs")
 scopes<-c("Overall disease groups","BMI25 four groups");groups<-c("Control","T2D","NC","OC","NT2D","OT2D")
 r1_need(setequal(unique(d$scope),scopes)&&setequal(unique(d$group),groups),"Composition scopes/groups differ")
 r1_need(all(is.finite(d$raw_abundance)&d$raw_abundance>=0),"Invalid saved abundance")
 sums<-aggregate(d$raw_abundance,list(group=d$group),sum);names(sums)[2]<-"saved_total_percent"
 r1_need(all(sums$saved_total_percent>if(species)93 else 99)&all(sums$saved_total_percent<if(species)97 else 100.1),"Saved totals differ; do not reclose")
 z<-d;z$scope<-r1_factor(z$scope,scopes);z$group<-r1_factor(z$group,groups)
 # Keep the saved category name verbatim. The palette's 'Other' key is a
 # style key, not the source data's phylum category ('Other phyla').
 # The current-run composition producer stores display_name as a factor.
 # c(factor_values, character_value) can concatenate internal integer codes;
 # form the display order from labels explicitly, leaving d unchanged.
 levels<-unique(as.character(d$display_name));other<-if(species)"Other retained features"else"Other phyla"
 r1_need(other%in%levels,"Saved Other category missing")
 levels<-c(levels[levels!=other],other);z$display_name<-r1_factor(z$display_name,levels)
 palette<-if(species)setNames(unlist(style$species_palette),levels)else unlist(style$phylum_palette)
 if(!species){
  r1_need(sum(names(palette)=="Other")==1L,"Phylum palette needs one explicit Other colour")
  names(palette)[names(palette)=="Other"]<-other
 }
 r1_need(length(levels)==if(species)9L else 5L,"Saved composition feature count differs")
 r1_need(all(levels%in%names(palette)),"Unmapped composition colour")
 label_fun<-function(x)as.expression(lapply(x,function(s)if(s==other)as.name(s)else call("italic",s)))
 p<-ggplot2::ggplot(z,ggplot2::aes(x=group,y=raw_abundance,fill=display_name))+
  ggplot2::geom_col(width=style$composition$bar_width,colour=style$composition$segment_border_colour,linewidth=style$composition$segment_border_linewidth)+
  ggplot2::facet_grid(cols=ggplot2::vars(scope),scales="free_x",space="free_x")+
  ggplot2::scale_fill_manual(values=palette,breaks=levels,labels=if(species)label_fun else ggplot2::waiver(),drop=FALSE)+
  ggplot2::scale_y_continuous(limits=c(0,100),expand=ggplot2::expansion(mult=c(0,.02)))+
  ggplot2::labs(x=NULL,y=style$composition$y_axis_title,fill=NULL)+theme_round1(style)+
  ggplot2::theme(legend.position=style$composition$legend_position,
   strip.background=ggplot2::element_rect(fill=style$composition$strip_fill,colour=NA))
 p<-r1_titles(p,panel,"",style)
 list(plot=p,display_data=d,audits=setNames(list(sums),paste0(panel$output_filename,"_saved_group_totals.csv")),
  caption=if(species)"Top 8 species and remaining retained features; values were not reclosed to 100%."else"Saved Top 4 phyla and Other; no abundance aggregation or renormalization was performed.")
}
