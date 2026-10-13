#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
root <- dirname(dirname(normalizePath(script)))
if (dir.exists(file.path(root, ".Rlib"))) .libPaths(c(file.path(root, ".Rlib"), .libPaths()))
source(file.path(root, "R", "receptor_nb.R"))
source(file.path(root, "R", "config.R"))
if (length(args) != 1L) stop("Usage: Rscript scripts/run_nb.R config.txt")
a <- read_analysis_config(args[1L])
if (dir.exists(a$output) || file.exists(a$output)) stop("Output already exists; choose a new directory")
config <- if (is.null(a$config)) list() else a$config
write_analysis <- function(result, output) {
  dir.create(output, recursive = TRUE)
  saveRDS(result, file.path(output, "analysis.rds"))
  results_file <- if (result$analysis_type == "gene_set_sum") "gene_set_results.csv" else "gene_results.csv"
  write.csv(result$results, file.path(output, results_file), row.names = FALSE)
  if (!is.null(result$gene_set_members)) {
    write.csv(result$gene_set_members, file.path(output, "gene_set_members.csv"), row.names = FALSE)
  }
  write.csv(result$effects, file.path(output, "receptor_effects.csv"), row.names = FALSE)
  writeLines(result$session, file.path(output, "sessionInfo.txt"))
  message("Wrote results to ", output, ": ", sum(result$results$status == "tested"), " tested, ",
          sum(result$results$status == "failed"), " failed, ", sum(result$results$status == "filtered"), " filtered")
  if (!any(result$results$status == "tested")) stop("No valid tests; inspect ", results_file, " for failure reasons")
}
if (is.null(a$receptors_file)) {
  result <- run_receptor_nb(a$input, a$receptors, a$cell_type, config)
  write_analysis(result, a$output)
} else {
  groups <- read_receptor_file(a$receptors_file)
  dir.create(a$output, recursive = TRUE)
  manifest <- data.frame(
    analysis_id = vapply(groups, function(g) sprintf("analysis_%04d", g$line), character(1)),
    line = vapply(groups, function(g) g$line, integer(1)),
    receptors = vapply(groups, function(g) paste(g$receptors, collapse = ","), character(1)),
    status = "pending", message = "")
  write.csv(manifest, file.path(a$output, "analyses.csv"), row.names = FALSE)
  for (i in seq_along(groups)) {
    message("Running ", manifest$analysis_id[i], ": ", manifest$receptors[i])
    error <- tryCatch({
      result <- run_receptor_nb(a$input, groups[[i]]$receptors, a$cell_type, config,
                                pool_receptors = TRUE)
      result$receptors_file <- normalizePath(a$receptors_file)
      result$receptors_line <- groups[[i]]$line
      write_analysis(result, file.path(a$output, manifest$analysis_id[i]))
      NULL
    }, error = conditionMessage)
    manifest$status[i] <- if (is.null(error)) "completed" else "failed"
    manifest$message[i] <- if (is.null(error)) "" else error
    write.csv(manifest, file.path(a$output, "analyses.csv"), row.names = FALSE)
  }
  if (any(manifest$status == "failed")) stop("One or more analyses failed; inspect analyses.csv")
}
