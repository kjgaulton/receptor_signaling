analysis <- list(
  input = "data.rds",
  output = "results/il7r_cd4_nb",
  receptors = c("IL7R"),
  cell_type = "CD4 T",
  config = list(
    cell_type_key = "cell_type",
    donor_key = "donor_id",
    assay = "RNA",
    layer = "counts",
    covariates = character(), # e.g. c("batch", "condition") if estimable
    mito_prefix = "MT-",      # mouse often "mt-"; alternatively use mito_genes
    family = "nbinom2",
    receptor_df = 1L,         # 3L for a joint nonlinear spline test
    qc_df = 3L,
    min_detected = 20L,
    min_fraction = 0.01,
    min_donors = 5L,
    min_receptor_donors = 3L,
    adjust_detected = TRUE
  )
)
