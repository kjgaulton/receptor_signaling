# Non-executable key=value configuration. Paths are relative to this file.
read_analysis_config <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  lines <- trimws(readLines(path, warn = FALSE))
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  values <- list()
  for (line in lines) {
    pos <- regexpr("=", line, fixed = TRUE)[1L]
    if (pos < 2L) stop("Expected key = value: ", line)
    key <- trimws(substr(line, 1L, pos - 1L))
    value <- trimws(substring(line, pos + 1L))
    if (key %in% names(values)) stop("Duplicate configuration key: ", key)
    values[[key]] <- value
  }
  top <- c("input", "output", "receptors", "receptors_file", "cell_type")
  allowed <- c(top, setdiff(names(nb_defaults()), "gene_sets"))
  unknown <- setdiff(names(values), allowed)
  if (length(unknown)) stop("Unknown configuration keys: ", paste(unknown, collapse = ", "))
  required <- c("input", "output", "cell_type")
  if (!all(required %in% names(values)) || any(!nzchar(unlist(values[required])))) {
    stop("Config requires nonempty input, output, cell_type")
  }
  if (sum(c("receptors", "receptors_file") %in% names(values)) != 1L) {
    stop("Specify exactly one of receptors_file or receptors")
  }
  if (any(!nzchar(unlist(values[intersect(c("receptors", "receptors_file"), names(values))])))) {
    stop("Receptor input must not be empty")
  }
  vectors <- c("receptors", "covariates", "mito_genes", "targets")
  numeric_keys <- c("receptor_df", "qc_df", "min_detected", "min_fraction", "min_donors", "min_receptor_donors")
  for (key in names(values)) {
    value <- values[[key]]
    if (key %in% vectors) {
      if (!nzchar(value)) {
        values[[key]] <- if (key == "covariates") character() else NULL
      } else {
        parts <- trimws(strsplit(paste0(value, ",__END__"), ",", fixed = TRUE)[[1L]])
        parts <- head(parts, -1L)
        if (any(!nzchar(parts))) stop("Empty list item: ", key)
        values[[key]] <- parts
      }
    } else if (key %in% numeric_keys) {
      number <- suppressWarnings(as.numeric(value))
      if (length(number) != 1L || !is.finite(number)) stop("Invalid numeric value: ", key)
      values[[key]] <- number
    } else if (key == "adjust_detected") {
      if (!tolower(value) %in% c("true", "false")) stop("adjust_detected must be true or false")
      values[[key]] <- tolower(value) == "true"
    } else if (!nzchar(value)) {
      if (key == "gene_sets_gmt") values[[key]] <- NULL else stop("Empty value: ", key)
    }
  }
  for (key in intersect(c("input", "output", "gene_sets_gmt", "receptors_file"), names(values))) {
    value <- path.expand(values[[key]])
    if (!grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", value)) value <- file.path(dirname(path), value)
    values[[key]] <- value
  }
  c(values[intersect(top, names(values))], list(config = values[setdiff(names(values), top)]))
}


# Preserve physical line numbers for traceable, filesystem-safe analysis IDs.
read_receptor_file <- function(path) {
  lines <- trimws(readLines(path, warn = FALSE))
  groups <- list()
  for (i in seq_along(lines)) {
    if (!nzchar(lines[i]) || startsWith(lines[i], "#")) next
    members <- trimws(head(strsplit(paste0(lines[i], ",__END__"), ",", fixed = TRUE)[[1L]], -1L))
    if (any(!nzchar(members))) stop("Empty gene identifier on receptor-file line ", i)
    groups[[length(groups) + 1L]] <- list(line = i, receptors = unique(members))
  }
  if (!length(groups)) stop("Receptor file contains no analyses")
  groups
}
