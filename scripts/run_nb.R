#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
root <- dirname(dirname(normalizePath(script)))
if (dir.exists(file.path(root, ".Rlib"))) .libPaths(c(file.path(root, ".Rlib"), .libPaths()))
source(file.path(root, "R", "receptor_nb.R"))
if (length(args) != 1L) stop("Usage: Rscript scripts/run_nb.R config.R")
env <- new.env(parent = baseenv())
sys.source(args[1L], envir = env)
if (!exists("analysis", envir = env, inherits = FALSE)) stop("Config must define an analysis list")
a <- env$analysis
if (!all(c("input", "receptors", "cell_type", "output") %in% names(a))) stop("Config requires input, receptors, cell_type, output")
if (dir.exists(a$output) || file.exists(a$output)) stop("Output already exists; choose a new directory")
config <- if (is.null(a$config)) list() else a$config
result <- run_receptor_nb(a$input, a$receptors, a$cell_type, config)
dir.create(a$output, recursive = TRUE)
saveRDS(result, file.path(a$output, "analysis.rds"))
write.csv(result$results, file.path(a$output, "gene_results.csv"), row.names = FALSE)
write.csv(result$effects, file.path(a$output, "receptor_effects.csv"), row.names = FALSE)
writeLines(result$session, file.path(a$output, "sessionInfo.txt"))
message("Wrote results to ", a$output, ": ", sum(result$results$status == "tested"), " tested, ",
        sum(result$results$status == "failed"), " failed, ", sum(result$results$status == "filtered"), " filtered")
if (!any(result$results$status == "tested")) stop("No valid tests; inspect gene_results.csv for failure reasons")
