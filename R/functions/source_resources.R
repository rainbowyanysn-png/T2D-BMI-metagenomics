sha256_file <- function(path) digest::digest(path, algo = "sha256", file = TRUE)

hub_cache_info <- function() {
  cache <- ExperimentHub::getExperimentHubOption("CACHE")
  database <- file.path(cache, "BiocFileCache.sqlite")
  if (!file.exists(database)) return(data.frame())
  connection <- DBI::dbConnect(RSQLite::SQLite(), dbname = database, flags = RSQLite::SQLITE_RO)
  on.exit(DBI::dbDisconnect(connection))
  DBI::dbReadTable(connection, "resource")
}

cached_resource_path <- function(id, info = hub_cache_info()) {
  if (!nrow(info)) return(NULL)
  hit <- info[info$rname == id | startsWith(info$rname, paste0(id, " :")), , drop = FALSE]
  if (nrow(hit) != 1L) return(NULL)
  p <- hit$rpath[[1]]
  if (!grepl("^([A-Za-z]:[/\\\\]|/|\\\\\\\\)", p)) p <- file.path(ExperimentHub::getExperimentHubOption("CACHE"), p)
  if (file.exists(p)) p else NULL
}

prepare_frozen_resources <- function(run, config) {
  previous_timeout <- getOption("timeout")
  options(timeout = max(previous_timeout, config$source_data$download_timeout_seconds))
  on.exit(options(timeout = previous_timeout), add = TRUE)
  require_project_packages(c("ExperimentHub", "S4Vectors", "DBI", "RSQLite", "digest"))
  object_manifest <- utils::read.csv(file.path(run$root, "config/package_input_objects.csv"), check.names = FALSE)
  for (i in seq_len(nrow(object_manifest))) {
    name <- object_manifest$object[i]
    object <- get0(name, asNamespace("curatedMetagenomicData"), inherits = TRUE)
    if (is.null(object)) object <- getExportedValue("curatedMetagenomicData", name)
    hash <- digest::digest(object, algo = "sha256", serialize = TRUE, serializeVersion = 2)
    assert_true(hash == object_manifest$sha256[i], paste("Packaged raw-input content mismatch:", name))
  }
  manifest_file <- file.path(run$root, config$source_data$resource_manifest)
  assert_file(manifest_file, "Frozen source-resource manifest")
  manifest <- utils::read.csv(manifest_file, check.names = FALSE)
  assert_true(nrow(manifest) > 0 && !anyDuplicated(manifest$resource_id), "Invalid resource manifest")
  local <- isTRUE(config$source_data$prefer_local_cache)
  inventory <- hub_cache_info()
  catalog_cached <- nrow(inventory) > 0L && any(inventory$rname == "experimenthub.sqlite3")
  if (!catalog_cached) {
    assert_true(isTRUE(config$source_data$allow_download_if_missing), "No cached ExperimentHub catalog and downloading is disabled.")
    local <- FALSE
  }
  cache_errors <- character()
  hub <- tryCatch(ExperimentHub::ExperimentHub(localHub = local), error = function(e) {
    cache_errors <<- c(cache_errors, paste("Local catalog:", conditionMessage(e)))
    if (!isTRUE(config$source_data$allow_download_if_missing)) stop(e)
    tryCatch(ExperimentHub::ExperimentHub(localHub = FALSE), error = function(network_error) {
      safe_write_lines(c(cache_errors, paste("Network catalog:", conditionMessage(network_error))), path_in_run(run, "logs", "resource_acquisition_failure.txt"))
      stop("Cannot obtain the frozen ExperimentHub resource catalog. See resource_acquisition_failure.txt; no replacement data have been used.", call. = FALSE)
    })
  })
  metadata <- as.data.frame(S4Vectors::mcols(hub))
  missing <- setdiff(manifest$resource_id, rownames(metadata))
  if (length(missing) && isTRUE(config$source_data$allow_download_if_missing)) {
    hub <- ExperimentHub::ExperimentHub(localHub = FALSE)
    metadata <- as.data.frame(S4Vectors::mcols(hub))
  }
  assert_true(all(manifest$resource_id %in% rownames(metadata)), "Frozen resource IDs unavailable. No substitute resources permitted.")
  assert_true(identical(as.character(metadata[manifest$resource_id, "title"]), manifest$title), "Frozen resource ID/title mismatch.")
  actual <- vector("list", nrow(manifest))
  for (i in seq_len(nrow(manifest))) {
    id <- manifest$resource_id[i]
    path <- cached_resource_path(id)
    downloaded <- FALSE
    if (is.null(path)) {
      assert_true(isTRUE(config$source_data$allow_download_if_missing), paste("Frozen resource missing from cache:", id))
      # Request this exact ID; never re-select a newer study/data-type resource.
      attempts <- list()
      acquired <- FALSE
      for (attempt in seq_len(config$source_data$download_max_attempts)) {
        message("Frozen resource download: ", id, " attempt ", attempt)
        acquired <- tryCatch({
          online_hub <- ExperimentHub::ExperimentHub(localHub = FALSE)
          temporary <- online_hub[[id]]
          rm(temporary)
          TRUE
        }, error = function(e) {
          attempts[[length(attempts) + 1L]] <<- data.frame(resource_id = id, attempt = attempt, time_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), error = conditionMessage(e))
          FALSE
        })
        if (acquired) break
      }
      if (length(attempts)) safe_write_csv(do.call(rbind, attempts), path_in_run(run, "logs", paste0("download_attempts_", id, ".csv")))
      assert_true(acquired, paste("Frozen resource unavailable after bounded retries:", id, "; no substitute version used. See download attempt log."))
      downloaded <- TRUE
      path <- cached_resource_path(id)
    }
    assert_true(!is.null(path), paste("Cannot locate exact resource after acquisition:", id))
    hash <- sha256_file(path)
    assert_true(tolower(hash) == tolower(manifest$sha256[i]), paste("SHA-256 mismatch for frozen raw resource:", id))
    actual[[i]] <- data.frame(resource_id = id, title = manifest$title[i], sha256 = hash, bytes = file.info(path)$size, downloaded_this_run = downloaded, checked_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
  }
  safe_write_csv(do.call(rbind, actual), path_in_run(run, "logs", "source_resource_verification.csv"))
  ExperimentHub::setExperimentHubOption("LOCAL", TRUE)
  pinned_hub <- ExperimentHub::ExperimentHub(localHub = TRUE)[manifest$resource_id]
  # Reuse the exact installed 3.16.1 processing code, injecting only a restricted Hub.
  # This local environment does not modify the installed package namespace.
  provider <- get("curatedMetagenomicData", asNamespace("curatedMetagenomicData"))
  retrieval <- get("returnSamples", asNamespace("curatedMetagenomicData"))
  scope <- new.env(parent = environment(provider))
  scope$ExperimentHub <- function(...) pinned_hub
  environment(provider) <- scope
  scope$curatedMetagenomicData <- provider
  environment(retrieval) <- scope
  retrieval
}
