publication_panel_theme <- function(config) {
  ggplot2::theme_classic(base_family = config$output$font_family, base_size = config$publication$font_size) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      axis.text = ggplot2::element_text(colour = "black"),
      axis.title = ggplot2::element_text(face = "bold"), legend.title = ggplot2::element_text(size = 8),
      legend.background=ggplot2::element_blank(),legend.box.background=ggplot2::element_blank(),legend.key=ggplot2::element_blank(),
      plot.caption = ggplot2::element_text(hjust = 0, size = 7),
      strip.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(6, 8, 6, 8))
}
