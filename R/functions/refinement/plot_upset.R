plot_upset <- function(d,panel,root,paths,style) {
 sets<-unlist(style$round2$comparisons);r2_columns(d,c("feature",sets,"intersection_pattern"))
 counts<-r1_read(root,paths,"Fig3C_intersection_counts.csv");sizes<-r1_read(root,paths,"Fig3C_set_sizes.csv")
 r2_columns(counts,c("intersection","n","pattern"));r2_columns(sizes,c("comparison","n"))
 r1_need(nrow(d)==79L&&!anyDuplicated(d$feature)&&sum(counts$n)==79L,"UpSet union differs")
 r1_need(identical(as.character(sizes$comparison),sets)&&identical(as.integer(sizes$n),c(69L,32L,10L,4L,5L)),"Saved set sizes differ")
 r1_need(!anyDuplicated(counts$pattern)&&!anyDuplicated(counts$intersection),"Duplicate saved intersection")
 r1_need(identical(as.integer(colSums(as.matrix(d[sets]))),as.integer(sizes$n)),"Membership/set-size inconsistency")
 observed<-table(factor(d$intersection_pattern,levels=counts$pattern))
 r1_need(!anyNA(match(d$intersection_pattern,counts$pattern))&&identical(as.integer(observed),as.integer(counts$n)),"Saved intersection counts disagree with memberships")
 bits<-strsplit(counts$pattern,"|",fixed=TRUE)
 r1_need(all(lengths(bits)==5L)&&all(unlist(bits)%in%c("TRUE","FALSE")),"Invalid saved membership patterns")
 display_sets<-unlist(style$upset$display_order)
 r1_need(identical(display_sets,c("OC vs NC","OT2D vs NT2D","OT2D vs OC","NT2D vs NC","T2D vs Control")),"UpSet display order differs")
 ypos<-setNames(5:1,display_sets)
 # Decode the saved patterns in ORIGINAL set order; remap y only for display.
 dots<-do.call(rbind,lapply(seq_along(bits),function(i)data.frame(intersection=counts$intersection[i],comparison=sets,y=unname(ypos[sets]),member=bits[[i]]=="TRUE")))
 segments<-do.call(rbind,lapply(counts$intersection,function(i){y<-dots$y[dots$intersection==i&dots$member];data.frame(intersection=i,lo=min(y),hi=max(y))}))
 s<-style$upset
 r1_need(all(unlist(s[c("set_panel_width_mm","name_panel_width_mm","matrix_height_mm")])>0)&&s$column_gap_mm>=0,"Invalid compact UpSet spacing")
 theme<-theme_round1(style)+ggplot2::theme(legend.position="none",plot.title=ggplot2::element_blank(),
  axis.text.x=ggplot2::element_text(margin=ggplot2::margin(t=1,unit="mm")),
  axis.title.x=ggplot2::element_text(margin=ggplot2::margin(t=1,unit="mm")),
  plot.margin=ggplot2::margin(0,0,0,0,unit="mm"))
 top<-ggplot2::ggplot(counts,ggplot2::aes(x=intersection,y=n))+
  ggplot2::geom_col(width=s$intersection_bar_width,fill=s$intersection_bar_colour)+
  ggplot2::geom_text(ggplot2::aes(label=n),vjust=-.4,family=style$font$family,size=2.8)+
  ggplot2::scale_x_continuous(breaks=counts$intersection,limits=range(counts$intersection)+c(-.5,.5),expand=c(0,0))+
  ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(0,.15)))+
  ggplot2::labs(x=NULL,y="Intersection size")+theme+
  ggplot2::theme(axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank(),axis.line.x=ggplot2::element_blank(),plot.margin=ggplot2::margin(0,0,s$row_margin_mm,0,unit="mm"))
 matrix<-ggplot2::ggplot(dots,ggplot2::aes(x=intersection,y=y))+
  ggplot2::geom_segment(data=segments,ggplot2::aes(x=intersection,xend=intersection,y=lo,yend=hi),inherit.aes=FALSE,colour=s$matrix_line_colour,linewidth=.4)+
  ggplot2::geom_point(ggplot2::aes(colour=member),size=s$member_point_size)+
  ggplot2::scale_colour_manual(values=c(`TRUE`=s$member_colour,`FALSE`=s$nonmember_colour))+
  ggplot2::scale_x_continuous(breaks=counts$intersection,limits=range(counts$intersection)+c(-.5,.5),expand=c(0,0))+
  ggplot2::scale_y_continuous(breaks=1:5,labels=rev(display_sets),limits=c(.5,5.5),expand=c(0,0))+
  ggplot2::labs(x="Intersection",y=NULL)+theme+
  ggplot2::theme(axis.text.y=ggplot2::element_blank(),axis.ticks.y=ggplot2::element_blank(),axis.line.y=ggplot2::element_blank(),plot.margin=ggplot2::margin(s$row_margin_mm,0,0,s$column_gap_mm,unit="mm"))
 sizes$y<-unname(ypos[sizes$comparison])
 left<-ggplot2::ggplot(sizes,ggplot2::aes(y=y,x=n))+
  ggplot2::geom_col(orientation="y",width=s$set_bar_width,fill=s$set_bar_colour)+
  ggplot2::geom_text(ggplot2::aes(label=n),hjust=1.15,family=style$font$family,size=2.7,colour=s$set_bar_colour)+
  ggplot2::scale_x_reverse(expand=ggplot2::expansion(mult=c(.08,.01)))+
  ggplot2::scale_y_continuous(breaks=1:5,labels=rev(display_sets),limits=c(.5,5.5),expand=c(0,0))+
  ggplot2::labs(x="Set size",y=NULL)+theme+
  ggplot2::theme(axis.line.y=ggplot2::element_blank(),axis.ticks.y=ggplot2::element_blank(),axis.text.y=ggplot2::element_blank(),plot.margin=ggplot2::margin(s$row_margin_mm,0,0,0,unit="mm"))
 # A separate aligned middle panel prevents names appearing outside the set bars.
 name_data<-data.frame(y=5:1,label=display_sets)
 names_panel<-ggplot2::ggplot(name_data,ggplot2::aes(x=1,y=y,label=label))+
  ggplot2::geom_text(family=style$font$family,size=s$name_size,hjust=1)+
  ggplot2::scale_x_continuous(limits=c(0,1),expand=c(0,0))+
  ggplot2::scale_y_continuous(limits=c(.5,5.5),expand=c(0,0))+
  theme+ggplot2::theme(axis.text=ggplot2::element_blank(),axis.ticks=ggplot2::element_blank(),axis.line=ggplot2::element_blank(),axis.title=ggplot2::element_blank(),plot.margin=ggplot2::margin(s$row_margin_mm,0,0,0,unit="mm"))
 # Fixed physical label width avoids the large empty gaps of a proportional column.
 # Fixed bar/name widths and zero proportional gutters put labels next to dots.
 widths<-grid::unit(c(s$set_panel_width_mm,s$name_panel_width_mm,1),c("mm","mm","null"))
 heights<-grid::unit(c(1,s$matrix_height_mm),c("null","mm"))
 p<-patchwork::wrap_plots(patchwork::plot_spacer(),patchwork::plot_spacer(),top,left,names_panel,matrix,ncol=3,widths=widths,heights=heights)+
  patchwork::plot_annotation(title=panel$title,theme=ggplot2::theme(text=ggplot2::element_text(family=style$font$family),plot.title=ggplot2::element_text(size=style$font$title_size,hjust=.5,face="bold"),plot.background=ggplot2::element_rect(fill="white",colour=NA)))
 list(plot=p,display_data=d,audits=list(Fig3C_display_spacing.csv=data.frame(set_panel_width_mm=s$set_panel_width_mm,name_panel_width_mm=s$name_panel_width_mm,column_gap_mm=s$column_gap_mm,matrix_height_mm=s$matrix_height_mm),Fig3C_saved_intersection_counts.csv=counts,Fig3C_saved_set_sizes.csv=sizes,
  Fig3C_display_order.csv=data.frame(position=seq_along(display_sets),comparison=display_sets,n=sizes$n[match(display_sets,sizes$comparison)])),
  caption="Memberships, intersection counts and set sizes are saved RC3 inputs. Union: 79 species. Original set sizes remain 69/32/10/4/5. Display rows (top to bottom): OC vs NC, OT2D vs NT2D, OT2D vs OC, NT2D vs NC, T2D vs Control. Group labels are between the set-size bars and intersection matrix. No significance set was regenerated.")
}
