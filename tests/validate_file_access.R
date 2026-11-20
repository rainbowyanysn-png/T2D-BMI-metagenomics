# R-level fail-fast interception complements static dependency inspection.
# This is not an operating-system filesystem sandbox.
install_legacy_access_guard <- function(prefixes) {
  assert_true(length(prefixes) > 0L && all(nzchar(prefixes)), "Supply explicit legacy deny prefixes for clean acceptance")
  normalize <- function(x) tolower(gsub("\\\\", "/", x))
  prefixes <- sub("/+$", "", normalize(prefixes))
  events <- new.env(parent = emptyenv()); events$rows <- list()
  guard <- function(path) {
    if (!is.character(path) || !length(path)) return(invisible(NULL))
    text <- normalize(path)
    denied <- vapply(text, function(p) any(p == prefixes | startsWith(p, paste0(prefixes, "/"))), logical(1))
    if (any(denied)) {
      events$rows[[length(events$rows) + 1L]] <- data.frame(path = path[denied], denied = TRUE)
      stop("Clean-run guard blocked legacy filesystem access", call. = FALSE)
    }
    invisible(NULL)
  }
  assign("v3_guard", guard, envir = .GlobalEnv)
  assign("v3_file_access_report", function() if (length(events$rows)) do.call(rbind, events$rows) else data.frame(path = "No denied legacy-path access attempted", denied = FALSE), envir = .GlobalEnv)
  for (name in c("file", "load", "readRDS", "source", "sys.source", "readLines")) {
    argument <- if (name == "file") "description" else if (name == "readLines") "con" else "file"
    expression <- substitute(get("v3_guard", envir = .GlobalEnv)(ARGUMENT), list(ARGUMENT = as.name(argument)))
    trace(name, tracer = expression, print = FALSE, where = baseenv())
  }
  invisible(TRUE)
}
