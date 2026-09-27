# Raw-count, within-donor receptor association tests.

read_receptor_rds <- function(path, assay = "RNA", layer = "counts") {
  object <- readRDS(path)
  if (inherits(object, "Seurat")) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("Install SeuratObject")
    a <- object[[assay]]
    if (is.null(a)) stop("Assay not found: ", assay)
    if (inherits(a, "Assay5")) {
      layers <- SeuratObject::Layers(a, search = NA)
      if (!layer %in% layers) {
        stop("Exact count layer not found. For split Seurat v5 layers, first use ",
             "SeuratObject::JoinLayers(object, assay=..., layers='counts', new='counts').")
      }
      counts <- SeuratObject::LayerData(a, layer = layer)
    } else {
      counts <- SeuratObject::GetAssayData(a, slot = layer)
    }
    metadata <- object[[]]
  } else if (inherits(object, "SingleCellExperiment") || inherits(object, "SummarizedExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) stop("Install SummarizedExperiment")
    counts <- SummarizedExperiment::assay(object, layer)
    metadata <- as.data.frame(SummarizedExperiment::colData(object))
  } else if (is.list(object) && all(c("counts", "metadata") %in% names(object))) {
    counts <- object$counts
    metadata <- as.data.frame(object$metadata)
  } else stop("RDS must contain Seurat, SingleCellExperiment, or list(counts=..., metadata=...)")
  if (!(is.matrix(counts) || inherits(counts, "Matrix"))) {
    stop("Counts must be an in-memory matrix or Matrix sparse matrix")
  }
  if (is.null(colnames(counts)) || is.null(rownames(counts)) ||
      anyDuplicated(colnames(counts)) || anyDuplicated(rownames(counts)) ||
      anyDuplicated(rownames(metadata)) || !all(colnames(counts) %in% rownames(metadata))) {
    stop("Unique gene/cell names and matching metadata row names are required")
  }
  list(counts = counts, metadata = metadata[colnames(counts), , drop = FALSE])
}

nb_defaults <- function() list(
  cell_type_key = "cell_type", donor_key = "donor_id", assay = "RNA", layer = "counts",
  covariates = character(), mito_prefix = "MT-", mito_genes = NULL,
  family = "nbinom2", receptor_df = 1L, qc_df = 3L,
  min_detected = 20L, min_fraction = 0.01, min_donors = 5L,
  min_receptor_donors = 3L, adjust_detected = TRUE, targets = NULL
)

scaled_basis <- function(x, df) {
  k <- length(unique(x))
  if (k < 2L || sd(x) < 1e-10) return(matrix(numeric(), length(x), 0L))
  z <- as.numeric(scale(x))
  if (df == 1L || k < df + 2L) return(matrix(z, ncol = 1L))
  splines::ns(z, df = df)
}

fit_nb_checked <- function(formula, data, family) {
  warnings <- character()
  fit <- withCallingHandlers(
    glmmTMB::glmmTMB(formula, data = data, family = family, REML = FALSE,
                    na.action = na.fail,
                    control = glmmTMB::glmmTMBControl(
                      optCtrl = list(iter.max = 1000L, eval.max = 1500L))),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  if (fit$fit$convergence != 0L || !isTRUE(fit$sdr$pdHess) || !is.finite(as.numeric(logLik(fit)))) {
    stop("Fit failed convergence / positive-definite Hessian checks; ", paste(warnings, collapse = "; "))
  }
  list(model = fit, warnings = warnings)
}

run_receptor_nb <- function(input, receptors, cell_type, config = list()) {
  if (!requireNamespace("glmmTMB", quietly = TRUE) || !requireNamespace("Matrix", quietly = TRUE)) {
    stop("Install glmmTMB and Matrix")
  }
  unknown <- setdiff(names(config), names(nb_defaults()))
  if (length(unknown)) stop("Unknown configuration fields: ", paste(unknown, collapse = ", "))
  cfg <- utils::modifyList(nb_defaults(), config)
  for (key in c("receptor_df", "qc_df", "min_detected", "min_donors", "min_receptor_donors")) {
    x <- cfg[[key]]
    if (length(x) != 1L || !is.finite(x) || x < 1 || x != as.integer(x)) stop("Invalid ", key)
  }
  if (cfg$min_donors < 3L) stop("At least three donors are required; five or more recommended")
  if (!is.finite(cfg$min_fraction) || cfg$min_fraction < 0 || cfg$min_fraction > 1) stop("Invalid min_fraction")
  if (!cfg$family %in% c("nbinom2", "nbinom1", "poisson")) stop("family must be nbinom2, nbinom1 or poisson")
  if (cfg$donor_key %in% cfg$covariates || cfg$cell_type_key %in% cfg$covariates) stop("Do not include donor/cell type as covariates")
  data <- read_receptor_rds(input, cfg$assay, cfg$layer)
  counts <- data$counts
  metadata <- data$metadata
  required <- c(cfg$cell_type_key, cfg$donor_key, cfg$covariates)
  if (!all(required %in% names(metadata))) stop("Missing metadata: ", paste(setdiff(required, names(metadata)), collapse = ", "))
  keep <- !is.na(metadata[[cfg$cell_type_key]]) & metadata[[cfg$cell_type_key]] == cell_type
  counts <- counts[, keep, drop = FALSE]
  metadata <- metadata[keep, , drop = FALSE]
  if (!ncol(counts)) stop("No cells match the requested cell type")
  if (anyNA(metadata[, required, drop = FALSE])) stop("Missing covariate/donor values; resolve before fitting")
  if (!inherits(counts, "sparseMatrix")) counts <- Matrix::Matrix(counts, sparse = TRUE)
  counts <- as(counts, "dgCMatrix")
  if (any(!is.finite(counts@x)) || any(counts@x < 0) || any(abs(counts@x - round(counts@x)) > 1e-8)) stop("Supply raw nonnegative integer counts")
  receptors <- unique(receptors)
  genes <- rownames(counts)
  if (!length(receptors) || !all(receptors %in% genes)) stop("Receptors must match unique count-matrix row names")
  donor <- droplevels(factor(metadata[[cfg$donor_key]]))
  if (nlevels(donor) < cfg$min_donors) stop("Insufficient donors for configured min_donors")
  mito <- if (is.null(cfg$mito_genes)) startsWith(genes, cfg$mito_prefix) else genes %in% cfg$mito_genes
  if (!any(mito)) stop("No mitochondrial genes matched; supply mito_prefix or mito_genes")
  ri <- match(receptors, genes)
  r <- t(as.matrix(counts[ri, , drop = FALSE]))
  for (k in seq_along(receptors)) {
    varying <- tapply(r[, k], donor, function(x) length(unique(x)) > 1L)
    if (sum(r[, k] > 0) < cfg$min_detected || sum(varying) < cfg$min_receptor_donors) {
      stop("Insufficient receptor detection/within-donor variation: ", receptors[k])
    }
  }
  total <- Matrix::colSums(counts) - rowSums(r)
  ndetected <- Matrix::colSums(counts > 0) - rowSums(r > 0)
  mt <- Matrix::colSums(counts[mito, , drop = FALSE]) - rowSums(r[, mito[ri], drop = FALSE])
  targets <- if (is.null(cfg$targets)) setdiff(genes, receptors) else unique(cfg$targets)
  if (!length(targets) || !all(targets %in% genes) || any(targets %in% receptors)) stop("Invalid target gene list")
  extra <- data.frame(row.names = seq_len(ncol(counts)))
  covariate_map <- character()
  for (k in seq_along(cfg$covariates)) {
    v <- metadata[[cfg$covariates[k]]]
    if (is.numeric(v)) {
      if (any(!is.finite(v))) stop("Nonfinite covariate: ", cfg$covariates[k])
      if (sd(v) < 1e-10) stop("Constant covariate: ", cfg$covariates[k])
      v <- as.numeric(scale(v))
    } else {
      v <- droplevels(factor(v))
      if (nlevels(v) < 2L) stop("Constant covariate: ", cfg$covariates[k])
    }
    name <- paste0("cov", k)
    extra[[name]] <- v
    covariate_map[name] <- cfg$covariates[k]
  }
  family <- switch(cfg$family, nbinom2 = glmmTMB::nbinom2(), nbinom1 = glmmTMB::nbinom1(), poisson = poisson())
  results <- vector("list", length(targets))
  effects <- list()
  for (i in seq_along(targets)) {
    gene <- targets[i]
    y <- as.numeric(counts[gene, ])
    row <- data.frame(gene = gene, status = "filtered", n_cells = length(y),
                      n_donors = nlevels(donor), n_detected = sum(y > 0),
                      lrt = NA_real_, df = NA_real_, p_value = NA_real_,
                      dispersion = NA_real_, donor_sd = NA_real_, warnings = "", message = "")
    if (sum(y > 0) < max(cfg$min_detected, cfg$min_fraction * length(y))) {
      row$message <- "Below detection threshold"
      results[[i]] <- row
      next
    }
    results[[i]] <- tryCatch({
      depth <- total - y
      if (any(depth <= 0)) stop("Nonpositive reference depth; supply whole-transcriptome counts and filter empty cells")
      d <- data.frame(y = y, donor = donor, log_depth = log(depth), extra)
      nuisance <- names(extra)
      qc <- list(mito = (mt - y * mito[match(gene, genes)]) / depth)
      if (cfg$adjust_detected) qc$detected <- log1p(ndetected - (y > 0))
      for (name in names(qc)) {
        b <- scaled_basis(qc[[name]], cfg$qc_df)
        for (j in seq_len(ncol(b))) {
          bn <- paste0("qc_", name, j)
          d[[bn]] <- b[, j]
          nuisance <- c(nuisance, bn)
        }
      }
      tested <- character()
      term_receptors <- character()
      for (k in seq_along(receptors)) {
        x <- log1p(1e4 * r[, k] / depth)
        # Keep the linear term on its natural log1p expression scale.
        b <- if (cfg$receptor_df == 1L) matrix(x, ncol = 1L) else scaled_basis(x, cfg$receptor_df)
        if (ncol(b) != cfg$receptor_df) stop("Insufficient expression diversity for requested receptor spline")
        for (j in seq_len(ncol(b))) {
          between <- ave(b[, j], donor, FUN = mean)
          within <- b[, j] - between
          if (sd(within) < 1e-10) stop("No within-donor receptor variation")
          wn <- paste0("r", k, "_within", j)
          bn <- paste0("r", k, "_between", j)
          d[[wn]] <- within
          tested <- c(tested, wn)
          term_receptors[wn] <- receptors[k]
          if (sd(between) > 1e-10) {
            d[[bn]] <- between
            nuisance <- c(nuisance, bn)
          }
        }
      }
      rhs0 <- paste(c("1", nuisance, "offset(log_depth)", "(1 | donor)"), collapse = " + ")
      rhs1 <- paste(c(rhs0, tested), collapse = " + ")
      design <- model.matrix(reformulate(c(nuisance, tested)), d)
      if (qr(design)$rank < ncol(design)) stop("Rank-deficient design; correlated receptors or confounded covariates")
      null <- fit_nb_checked(as.formula(paste("y ~", rhs0)), d, family)
      full <- fit_nb_checked(as.formula(paste("y ~", rhs1)), d, family)
      ll0 <- logLik(null$model)
      ll1 <- logLik(full$model)
      statistic <- 2 * as.numeric(ll1 - ll0)
      df <- attr(ll1, "df") - attr(ll0, "df")
      if (!is.finite(statistic) || statistic < -1e-6 || df != length(tested)) stop("Invalid nested model comparison")
      row$status <- "tested"
      row$lrt <- max(0, statistic)
      row$df <- df
      row$p_value <- pchisq(row$lrt, df = df, lower.tail = FALSE)
      row$dispersion <- sigma(full$model)
      row$donor_sd <- unname(attr(glmmTMB::VarCorr(full$model)$cond$donor, "stddev")[1L])
      row$warnings <- paste(unique(c(null$warnings, full$warnings)), collapse = "; ")
      co <- summary(full$model)$coefficients$cond[tested, , drop = FALSE]
      effects[[length(effects) + 1L]] <- data.frame(
        gene = gene, receptor = unname(term_receptors[tested]), term = tested,
        estimate = co[, 1L], se = co[, 2L],
        wald_lower = co[, 1L] - qnorm(.975) * co[, 2L],
        wald_upper = co[, 1L] + qnorm(.975) * co[, 2L],
        rate_ratio = if (cfg$receptor_df == 1L) exp(co[, 1L]) else NA_real_)
      row
    }, error = function(e) {
      row$status <- "failed"
      row$message <- conditionMessage(e)
      row
    })
    if (i %% 100L == 0L) message("Processed ", i, " / ", length(targets), " targets")
  }
  results <- do.call(rbind, results)
  results$q_value <- NA_real_
  tested <- results$status == "tested"
  # Failed eligible fits count toward the testing family, without invented p-values.
  if (any(tested)) results$q_value[tested] <- p.adjust(results$p_value[tested], "BH", n = sum(results$status != "filtered"))
  results <- results[order(results$q_value, results$p_value, na.last = TRUE), ]
  list(results = results, effects = if (length(effects)) do.call(rbind, effects) else data.frame(),
       config = cfg, receptors = receptors, cell_type = cell_type,
       donor_cell_counts = table(donor), covariate_map = covariate_map,
       input = normalizePath(input), session = capture.output(sessionInfo()),
       interpretation = "Asymptotic within-donor joint receptor LRT; BH within this run; model assumptions required")
}
