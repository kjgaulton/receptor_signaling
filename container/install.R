options(repos = c(CRAN = Sys.getenv("CRAN_SNAPSHOT")), timeout = 1200, Ncpus = 2L)
# Compile TMB and glmmTMB together from the same snapshot, avoiding binary ABI drift.
install.packages(c("glmmTMB", "SeuratObject"), type = "source")
for (p in c("glmmTMB", "SeuratObject")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("CRAN installation failed: ", p)
}
# Select the release repositories explicitly; do not let current-date release
# detection switch this older R image to a different Bioconductor release.
bioc <- paste0("https://bioconductor.org/packages/", Sys.getenv("BIOCONDUCTOR_VERSION"))
options(repos = c(BioCsoft = paste0(bioc, "/bioc"),
                 BioCann = paste0(bioc, "/data/annotation"),
                 BioCexp = paste0(bioc, "/data/experiment"),
                 CRAN = Sys.getenv("CRAN_SNAPSHOT")))
install.packages(c("BiocVersion", "SingleCellExperiment", "SummarizedExperiment"), type = "source")
required <- c("Matrix", "TMB", "glmmTMB", "SeuratObject", "SingleCellExperiment", "SummarizedExperiment")
for (p in required) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Required package failed to install: ", p)
}
dir.create("/opt/receptor-signaling", recursive = TRUE, showWarnings = FALSE)
packages <- as.data.frame(installed.packages()[, c("Package", "Version", "Built")])
write.csv(packages, "/opt/receptor-signaling/packages.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "/opt/receptor-signaling/build-session.txt")
