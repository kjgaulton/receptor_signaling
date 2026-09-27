.libPaths(c(normalizePath(".Rlib"), .libPaths()))
source("R/receptor_nb.R")
set.seed(221)
n <- 1200L
donor <- rep(seq_len(12L), each = 100L)
capture <- exp(rnorm(n, 0, .3))
background <- matrix(rpois(60L * n, rep(capture * 30, each = 60L)), 60L, n)
depth <- colSums(background)
r <- rnbinom(n, mu = capture * exp(rnorm(n, 1, .65)), size = 5)
x <- log1p(1e4 * r / depth)
w <- x - ave(x, donor)
intercept <- rnorm(12L, 0, .5)[donor]
signal <- rnbinom(n, mu = depth / 200 * exp(.65 * w + intercept), size = 6)
null <- rnbinom(n, mu = depth / 200 * exp(intercept), size = 6)
counts <- rbind(R = r, signal = signal, null = null, zero = 0, background)
rownames(counts)[5:nrow(counts)] <- c("MT-A", paste0("B", 1:59))
colnames(counts) <- paste0("c", seq_len(n))
counts <- Matrix::Matrix(counts, sparse = TRUE)
metadata <- data.frame(donor_id = as.character(donor), cell_type = "T", row.names = colnames(counts))
path <- tempfile(fileext = ".rds")
saveRDS(list(counts = counts, metadata = metadata[n:1, ]), path)
loaded <- read_receptor_rds(path)
stopifnot(identical(rownames(loaded$metadata), colnames(counts)))
cfg <- list(targets = c("signal", "null", "zero"), qc_df = 1L)
res <- run_receptor_nb(path, "R", "T", cfg)
print(res$results)
sig <- subset(res$results, gene == "signal")
stopifnot(sig$status == "tested", sig$p_value < 1e-6, sig$q_value < .05,
          subset(res$results, gene == "zero")$status == "filtered",
          subset(res$effects, gene == "signal")$estimate > 0,
          all(res$results$q_value >= res$results$p_value, na.rm = TRUE))
nonlinear <- run_receptor_nb(path, "R", "T", list(targets = "signal", receptor_df = 3L, qc_df = 1L))
stopifnot(nonlinear$results$status == "tested", nonlinear$results$df == 3L,
          all(is.na(nonlinear$effects$rate_ratio)))
# Independent receptors produce a joint two-degree-of-freedom test.
two <- rbind(counts, R2 = rnbinom(n, mu = capture * 4, size = 3))
saveRDS(list(counts = two, metadata = metadata), path)
joint <- run_receptor_nb(path, c("R", "R2"), "T", list(targets = "signal", qc_df = 1L))
stopifnot(joint$results$status == "tested", joint$results$df == 2L)
# Normalized/noninteger matrices must not be silently fitted as raw counts.
invalid <- counts
invalid[1, 1] <- .5
saveRDS(list(counts = invalid, metadata = metadata), path)
error <- tryCatch(run_receptor_nb(path, "R", "T", cfg), error = conditionMessage)
stopifnot(is.character(error), grepl("raw nonnegative integer", error))
# Identical receptors are not separately estimable; no p-values should be reported.
duplicate <- rbind(counts, Rcopy = counts["R", ])
saveRDS(list(counts = duplicate, metadata = metadata), path)
bad <- run_receptor_nb(path, c("R", "Rcopy"), "T", list(targets = "signal", qc_df = 1L))
stopifnot(bad$results$status == "failed", is.na(bad$results$p_value))
# Actual RDS roundtrips for both common single-cell object formats.
if (requireNamespace("SeuratObject", quietly = TRUE)) {
  obj <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = metadata)
  saveRDS(obj, path)
  got <- read_receptor_rds(path)
  stopifnot(isTRUE(all.equal(as.matrix(got$counts), as.matrix(counts))))
}
if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  obj <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts), colData = metadata)
  saveRDS(obj, path)
  got <- read_receptor_rds(path)
  stopifnot(isTRUE(all.equal(as.matrix(got$counts), as.matrix(counts))))
}
# CLI exercise, including saved significance and metadata files.
saveRDS(list(counts = counts, metadata = metadata), path)
out <- tempfile()
config <- tempfile(fileext = ".R")
analysis <- list(input = path, output = out, receptors = "R", cell_type = "T", config = cfg)
writeLines(c("analysis <-", capture.output(dput(analysis))), config)
status <- system2(file.path(R.home("bin"), "Rscript"), c("scripts/run_nb.R", shQuote(config)))
stopifnot(status == 0L, all(file.exists(file.path(out, c("gene_results.csv", "receptor_effects.csv", "analysis.rds", "sessionInfo.txt")))))
unlink(c(path, config))
unlink(out, recursive = TRUE)
cat("All NB tests passed: signal, spline LRT, BH, rank failure, RDS readers, metadata alignment and CLI.\n")
