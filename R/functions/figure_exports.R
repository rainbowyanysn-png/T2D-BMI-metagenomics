read_presentation_config <- function(root) yaml::read_yaml(file.path(root, "config/presentation.yml"))

check_export_package_versions <- function(root, packages) {
  versions <- utils::read.csv(file.path(root, "config/export_package_versions.csv"), check.names = FALSE)
  packages <- unique(c(packages, versions$package))
  require_project_packages(packages)
  for (package in packages) {
    row <- versions[versions$package == package, , drop = FALSE]
    assert_equal(nrow(row), 1L, paste("Frozen export package", package))
    assert_equal(as.character(utils::packageDescription(package, fields = "Version")), row$version, paste("Export package version", package))
  }
}

publication_theme <- function(config) {
  ggplot2::theme_bw(base_size = config$output$font_size, base_family = config$output$font_family) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black", linewidth = config$output$line_width),
      axis.text = ggplot2::element_text(colour = "black"), axis.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_blank(), plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.title = ggplot2::element_text(face = "bold", hjust = .5),
      plot.tag = ggplot2::element_text(face = "bold", size = 16), strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "bold"))
}

species_plot_label <- function(x) {
  label <- gsub("[._]+", " ", sub("^species[.]", "", x))
  label <- gsub(" sp ", " sp. ", label, fixed = TRUE)
  label <- gsub(" CAG ", " CAG-", label, fixed = TRUE)
  paste0("italic('", gsub("'", "", label), "')")
}

save_publication_figure <- function(plot, name, run, config, width, height, caption) {
  prefix <- path_in_run(run, "figures", name)
  paths <- paste0(prefix, c(".pdf", ".tiff", ".png"))
  assert_true(!any(file.exists(paths)), "Refusing to overwrite a figure")
  ggplot2::ggsave(paths[1], plot, width = width, height = height, units = "in", bg = "white", device = grDevices::cairo_pdf)
  ggplot2::ggsave(paths[2], plot, width = width, height = height, units = "in", dpi = config$output$figure_dpi, bg = "white", device = grDevices::tiff, type = "cairo", compression = config$output$tiff_compression)
  ggplot2::ggsave(paths[3], plot, width = width, height = height, units = "in", dpi = config$output$figure_dpi, bg = "white", device = grDevices::png, type = "cairo")
  safe_write_lines(caption, paste0(prefix, "_caption.txt"))
  data.frame(figure_id = name, width_inches = width, height_inches = height, dpi = config$output$figure_dpi, pdf = paste0(name, ".pdf"), tiff = paste0(name, ".tiff"), png = paste0(name, ".png"))
}

current_result_path <- function(run, relative_path) {
  assert_true(!grepl("(^|[/\\\\])[.][.]([/\\\\]|$)", relative_path), "Result path must stay within this run")
  directories <- if (!is.null(run$result_dirs)) unname(run$result_dirs) else run$dir
  candidates <- file.path(directories, relative_path)
  hits <- unique(candidates[file.exists(candidates)])
  assert_true(length(hits) == 1L, paste("Expected exactly one active result in this run:", relative_path))
  hits[1]
}

read_current_result <- function(run, relative_path) utils::read.csv(current_result_path(run, relative_path), check.names = FALSE)
