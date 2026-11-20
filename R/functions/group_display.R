analysis_group_levels <- c("Non_Overweight_Control", "Overweight_Control", "Non_Overweight_T2D", "Overweight_T2D")
group_display <- c(Non_Overweight_Control = "NC", Overweight_Control = "OC", Non_Overweight_T2D = "NT2D", Overweight_T2D = "OT2D")
display_group <- function(x) {
  if (!all(as.character(x) %in% names(group_display))) stop("Unknown analysis group in display mapping.", call. = FALSE)
  unname(group_display[as.character(x)])
}
