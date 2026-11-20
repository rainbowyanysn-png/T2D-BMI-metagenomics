assert_true <- function(condition, message) if (!isTRUE(condition)) stop(message, call. = FALSE)
assert_equal <- function(observed, expected, label) assert_true(identical(observed, expected), paste0(label, ": expected ", expected, ", observed ", observed))
assert_columns <- function(x, columns, label) { absent <- setdiff(columns, names(x)); assert_true(!length(absent), paste0(label, " lacks: ", paste(absent, collapse = ", "))) }
assert_file <- function(path, label = "file") assert_true(file.exists(path), paste0(label, " does not exist: ", path))
safe_write_csv <- function(x, path) { assert_true(!file.exists(path), paste0("Refusing to overwrite: ", path)); dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE); utils::write.csv(x, path, row.names = FALSE, na = "") }
safe_write_lines <- function(x, path) { assert_true(!file.exists(path), paste0("Refusing to overwrite: ", path)); dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE); writeLines(x, path) }
