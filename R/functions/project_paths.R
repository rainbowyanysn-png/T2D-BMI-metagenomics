activate_existing_project_library <- function(root) {
  if(as.character(getRversion())!="4.5.1")stop("This frozen environment requires R 4.5.1",call.=FALSE)
  explicit <- Sys.getenv("T2DBMI_R_LIBRARY")
  if (nzchar(explicit)) {
    if (!dir.exists(explicit)) stop("T2DBMI_R_LIBRARY does not exist; restore packages before running.", call. = FALSE)
    .libPaths(c(explicit, .libPaths()))
  } else if (requireNamespace("renv", quietly = TRUE)) {
    library <- renv::paths$library(project = root)
    if (dir.exists(library)) .libPaths(c(library, .libPaths()))
  }
  invisible(.libPaths())
}

project_root <- function(script = NULL) {
  if (is.null(script)) {
    argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(argument) != 1L) stop("Use Rscript with an explicit entry script.", call. = FALSE)
    script <- sub("^--file=", "", argument)
  }
  if (!grepl("^([A-Za-z]:[/\\\\]|/|\\\\\\\\)", script)) script <- file.path(getwd(), script)
  root <- normalizePath(dirname(script),winslash="/",mustWork=TRUE)
  required <- c("config", "R", "tests")
  if (!all(dir.exists(file.path(root, required)))) stop("Entry script is outside the project root; required directories are absent.", call. = FALSE)
  root
}

run_context <- function(root = project_root()) {
  output_root <- file.path(root, "output")
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  id <- paste0("Run_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  dir <- file.path(output_root, id)
  if (dir.exists(dir)) stop("Refusing to overwrite output directory: ", dir, call. = FALSE)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  writeLines("RC4_CURRENT_PROJECT_RUN",file.path(dir,"RC4_RUN_IDENTITY.txt"))
  list(root = root, id = id, dir = dir,
       derived_dir = file.path(root, "data_derived", id))
}

path_in_run <- function(run, ...) {
  path <- file.path(run$dir, ...)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}

# Identity is engineering provenance, not a statistical parameter.
assert_rc4_run_identity <- function(run) {
  marker<-file.path(run$dir,"RC4_RUN_IDENTITY.txt")
  if(!file.exists(marker)||!identical(readLines(marker,warn=FALSE),"RC4_CURRENT_PROJECT_RUN"))
    stop("Only a Run created by RC4 may be resumed or exported; older project Runs are rejected.",call.=FALSE)
  invisible(TRUE)
}
