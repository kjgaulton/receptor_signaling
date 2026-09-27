# Receptor–gene association from single-cell RDS data

This repository provides code to test for association between the expression of a signaling receptor and downstream targets using single cell profiles.  

## Input and quick start

An RDS file is a container, not a particular data format. Supported contents:

- Seurat objects, using the raw `RNA` assay `counts` layer by default.
- SingleCellExperiment / SummarizedExperiment objects with a `counts` assay.
- `list(counts = counts, metadata = metadata)`, with a **genes × cells** matrix
  and metadata row names matching cell names. Sparse Matrix objects are supported.

Supply whole-transcriptome raw UMI counts, not log-normalized, scaled, integrated,
SCT-corrected or highly-variable-gene-only data. Gene names must be unique. Receptor
names must match row names exactly; no automatic symbol/Ensembl conversion occurs.
Cells should already pass appropriate quality/doublet filtering.

```r
# One-time dependencies (already installed locally for this project):
dir.create(".Rlib", showWarnings = FALSE)
install.packages(c("Matrix", "glmmTMB"), lib = ".Rlib",
                 repos = "https://cloud.r-project.org")
# SeuratObject is needed for Seurat RDS input; SingleCellExperiment and
# SummarizedExperiment are needed for those object classes.
```

Edit [config.example.R](config.example.R), then run from this directory:

```bash
Rscript scripts/run_nb.R config.example.R
```

Specify input/output paths, receptor(s), cell type and donor metadata column.
With multiple receptors, one joint hypothesis is tested per target gene.
For Seurat v5 objects with split count layers, explicitly join the layers first:

```r
object <- SeuratObject::JoinLayers(object, assay = "RNA",
                                  layers = "counts", new = "counts")
saveRDS(object, "joined.rds")
```

The reader requires an exact layer match; it does not silently select one sample.
It reorders metadata by cell name and checks for missing/duplicate identifiers.
On-disk matrix backends are not supported by this initial implementation.

## Statistical model

For each target gene g and cell i from donor d:

```
Y_gi ~ NegativeBinomial(mean = mu_gi, variance = mu_gi + mu_gi^2 / theta_g)

log(mu_gi) = log(reference_UMIs_i) + intercept_g
          + QC_terms_i + supplied_covariates_i
          + donor_mean_receptor_terms_d
          + beta_g * within_donor_receptor_i
          + donor_random_intercept_gd
```

The response is the **raw target count**. Log reference depth is an offset with
coefficient fixed to one, rather than merely another predictor. The default
`nbinom2` family permits quadratic mean–variance overdispersion. `nbinom1` and
`poisson` can be explicitly chosen in the configuration. There is no automatic
per-gene distribution selection, and zero inflation is not assumed merely because
single-cell counts contain many zeros. Family choice needs dataset diagnostics.

For each target, reference UMI depth, mitochondrial fraction and detected-gene
count exclude the target and every input receptor. This avoids directly using
the response/receptor counts inside those nuisance quantities; it does not remove
all compositional bias or uncertainty in estimated exposure.

Receptor expression is `log1p(10000 * receptor_count / reference_UMIs)`.
The linear default decomposes it into its donor mean and its deviation from that
mean. Both models include the donor mean; only the full model includes the
within-donor deviation. Thus the primary test asks:

> Within this cell type and donor, is receptor RNA expression associated with
> target expression after the specified adjustment?

This differs from predicting expression in unseen donors and does not test a
between-donor receptor effect. A donor random intercept addresses donor-level
clustering under the model assumptions; it does not make arbitrary cell-level
dependence harmless. Use biological individual IDs across repeated samples.
Repeated sample effects beyond donor and donor-specific receptor slopes are not
yet modeled. Strong slope heterogeneity requires an extended model.

QC terms use natural splines (default three degrees of freedom) for mitochondrial
fraction and log detected-gene count. Constant QC terms are omitted. Depth enters
through its offset, not an additional estimated depth spline. Extra metadata
covariates are numeric linear terms or categorical fixed effects. Rank deficiency
is reported as a failed fit, rather than silently dropping confounded terms.

## Formal testing and effect sizes

Fit both models by maximum likelihood with the same family, cells, QC terms and
random-effect structure. The likelihood-ratio statistic is twice their log
likelihood difference. Its asymptotic chi-squared reference distribution has
degrees of freedom equal to the number of added receptor coefficients.

- `receptor_df = 1L`: one coefficient per receptor; a linear effect on the log mean
  versus transformed receptor expression. The exponentiated coefficient is a
  conditional rate ratio per one-unit increase in within-donor log1p normalized
  receptor expression. It is **not** a rate ratio per raw UMI or automatically per
  doubling of expression.
- `receptor_df = 3L`: three natural-spline basis terms per receptor, allowing a
  nonlinear curve. Each basis is decomposed into donor mean and within-donor
  deviation, and all within-donor terms are tested jointly. Basis coefficients
  are output, but no single rate ratio or direction is assigned to the curve.
- Multiple receptors: a joint omnibus test conditional on all their donor-mean
  terms. This is not a test of each receptor independently or of receptor synergy.
  Individual conditional LRTs and receptor interaction terms are not implemented.

Convergence failures, non-positive-definite Hessians, rank deficiency and invalid
nested comparisons yield **no p-value**. Other optimizer warnings are retained.
Near-zero random-effect variance and distributional misspecification can still
make otherwise numerically successful fits unreliable; inspect diagnostics.

BH correction applies across eligible target genes **within this run**. Failed
eligible fits count in the correction denominator. Low-detection genes are marked
filtered and excluded. Filtering should be prespecified; this procedure does not
prove that detection filtering preserves null calibration. Multiple receptor sets
or cell types require a deliberately defined larger testing family; separately
computed q-values do not control a combined family.

P-values are formal **model-based, asymptotic** tests, not distribution-free
guarantees. BH relies on suitable null p-values and dependence conditions. With
few donors, sparse receptors, misspecification or heterogeneous slopes, inspect
simulation-based diagnostics and consider a null parametric bootstrap or an
appropriate alternative model. These calibration procedures are not yet automated.
The default requires five donors and raw receptor variation in at least three;
these are guardrails, not guarantees of adequate power or calibration.

Do not interpret association as receptor activation or causal signaling. Receptor
mRNA is noisy, and unmeasured cell states can explain coexpression.

## Configuration and outputs

- `covariates`: e.g. batch or condition; only include scientifically justified,
  estimable terms. Do not pass donor ID, receptor/target measurements, or existing
  total-count QC columns. QC features are recomputed internally.
- `mito_prefix`: `MT-` by default; mouse annotations often use `mt-`.
  Alternatively pass a `mito_genes` vector, including Ensembl IDs if needed.
- `adjust_detected = FALSE`: sensitivity analysis omitting detected-gene count,
  which can capture real biology as well as technical variation.
- `min_detected`, `min_fraction`: prespecified target detection thresholds.
- `targets`: optional target gene vector for a pilot, retaining whole-transcriptome
  counts for normalization. Omitting it tests all genes except input receptors.

Outputs:

| File | Contents |
|---|---|
| `gene_results.csv` | LRT, degrees of freedom, p/q-values, fit status, dispersion, donor SD, warnings |
| `receptor_effects.csv` | Within-donor coefficients, standard errors, unadjusted 95% Wald intervals; linear rate ratios |
| `analysis.rds` | Results, settings, receptor names, donor cell counts, covariate mapping and session info |
| `sessionInfo.txt` | R and package versions |

Genome-wide fitting is sequential and requires two mixed models per eligible
gene. The count matrix remains sparse, but selected cells and metadata are held
in memory. There is no checkpoint/resume or parallel fitting yet. Existing output
directories are not overwritten. If every fit fails, outputs are saved for
diagnosis and the command returns an error.

## Validation

```bash
Rscript tests/test_nb.R
```

Synthetic tests exercise known-signal recovery, nonlinear LRT degrees of freedom,
BH behavior, rank-deficient receptor rejection, metadata alignment, supported RDS
readers and the command-line output path. They are implementation checks, **not a
genome-wide false-positive calibration study**. No user dataset has been fitted.

## References

- [glmmTMB families and NB variance functions](https://glmmtmb.github.io/glmmTMB/reference/nbinom2.html)
- [glmmTMB fitting, offsets and random effects](https://glmmtmb.github.io/glmmTMB/reference/glmmTMB.html)
- [glmmTMB troubleshooting and inference caveats](https://glmmtmb.github.io/glmmTMB/articles/troubleshooting.html)
- [Seurat v5 layers](https://satijalab.org/seurat/articles/seurat5_essential_commands)
- [Zimmerman et al., single-cell pseudoreplication](https://www.nature.com/articles/s41467-021-21038-1)
