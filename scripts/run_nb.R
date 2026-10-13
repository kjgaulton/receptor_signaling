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
result <- run_receptor_nb(a$input, a$receptors, a$cell_type, config)
dir.create(a$output, recursive = TRUE)
saveRDS(result, file.path(a$output, "analysis.rds"))
results_file <- if (result$analysis_type == "gene_set_sum") "gene_set_results.csv" else "gene_results.csv"
write.csv(result$results, file.path(a$output, results_file), row.names = FALSE)
if (!is.null(result$gene_set_members)) {
  write.csv(result$gene_set_members, file.path(a$output, "gene_set_members.csv"), row.names = FALSE)
}
write.csv(result$effects, file.path(a$output, "receptor_effects.csv"), row.names = FALSE)
writeLines(result$session, file.path(a$output, "sessionInfo.txt"))
message("Wrote results to ", a$output, ": ", sum(result$results$status == "tested"), " tested, ",
        sum(result$results$status == "failed"), " failed, ", sum(result$results$status == "filtered"), " filtered")
if (!any(result$results$status == "tested")) stop("No valid tests; inspect ", results_file, " for failure reasons")
