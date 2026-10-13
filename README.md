# Receptor–gene association from single-cell RDS data

This repository provides code to test for association between the expression of a signaling receptor and downstream targets using single cell profiles.  

## Portable container

The R workflow is packaged in a Docker image with its R dependencies, Seurat/SCE
readers, and a built-in synthetic self-test. See [CONTAINER.md](CONTAINER.md) for
build/run commands on macOS, Windows, Linux and HPC, and [container/config.txt](container/config.txt)
for a configuration using mounted input/output paths.

```bash
docker build -t receptor-signaling:0.1.0 .
docker run --rm receptor-signaling:0.1.0 --self-test
```

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

Edit [config.example.txt](config.example.txt), then run from this directory:

```bash
Rscript scripts/run_nb.R config.example.txt
```

Specify input/output paths, receptor(s), cell type and donor metadata column.
Use `receptors_file` to run one analysis per line, pooling genes on the same line.
For Seurat v5 objects with split count layers, explicitly join the layers first:

```r
object <- SeuratObject::JoinLayers(object, assay = "RNA",
                                  layers = "counts", new = "counts")
saveRDS(object, "joined.rds")
```

The reader requires an exact layer match; it does not silently select one sample.
It reorders metadata by cell name and checks for missing/duplicate identifiers.
On-disk matrix backends are not supported by this initial implementation.

## Plain-text configuration

Configurations are data, not executable R scripts. See [config.example.txt](config.example.txt).
Use one `key = value` per line, without quotes. Blank lines and lines beginning
with `#` are ignored; inline comments are not supported. Unknown or duplicate
keys are errors. `input`, `output`, `receptors_file`, and `cell_type` are required for batch runs.
Other keys use the defaults described below. Lists (`receptors`, `covariates`,
`targets`, `mito_genes`) use commas; booleans use `true` or `false` and numbers
use ordinary numeric notation (no R `L` suffix). An empty `covariates` means none;
empty `targets`, `mito_genes`, or `gene_sets_gmt` select their defaults.
Relative input, output, receptor-file, and GMT paths resolve from the configuration file's
location, independent of the shell's working directory. Paths and cell-type
names can contain spaces. Existing R-script configs must be converted to this
format; the CLI no longer sources them. The R function API still accepts lists.

## Receptor input file and batch analyses

Set `receptors_file = receptors.txt` (inside a container, for example
`receptors_file = /input/receptors.txt`). The file contains one analysis per line:

```text
IL7R
IL2RA,IL2RB,IL2RG
```

The first analysis uses IL7R expression. The second sums the raw counts of IL2RA,
IL2RB and IL2RG per cell, then applies the usual log1p normalization and
within-donor decomposition to that **single pooled predictor**. This gives one
linear coefficient (or `receptor_df` spline coefficients) regardless of the number
of pooled genes. It does not estimate separate gene effects or receptor synergy.
High-expression genes can dominate the pool. Detection/variation requirements
apply to pooled counts, while normalization and QC exclude all individual pool
members. The genes need not be annotated receptors; identifiers must match the
count matrix. Target sets containing any pool member are rejected as before.

Whitespace is trimmed, blank lines and full-line `#` comments are skipped,
and duplicate genes within a line are counted once. Empty comma-separated items
are errors. Repeated lines remain separate analyses. See [receptors.example.txt](receptors.example.txt).

Outputs are written to `output/analysis_0001/`, etc., using the original file
line number. Each folder contains the usual results and `analysis.rds` with pool
membership, source file and line. `output/analyses.csv` maps these folders to their
receptors and completion/failure messages. Failed analyses do not stop later
lines; the command exits nonzero if any line fails. Individual failed target fits
remain recorded in their result tables. Analyses run sequentially and reload the
input per line. The output parent must not already exist.

BH adjustment is performed separately within each line's analysis, not across
all receptor groups. It does not control a combined batch-wide testing family.
The same targets or GMT file and cell type are used for every line.

Legacy `receptors = IL7R` remains supported for a single analysis in the output
root. A legacy comma-separated `receptors` value retains the previous joint test
with separate predictors. Do not specify both `receptors` and `receptors_file`.
The R API supports pooling through `run_receptor_nb(..., pool_receptors = TRUE)`.

## Receptor–gene-set association

To test a pathway or an unsigned GRN target set across cells of one cell type,
set `gene_sets_gmt` in your plain-text configuration:

```text
gene_sets_gmt = pathways.gmt
```

The GMT file contains one set per line: set name, description (or `na`), then
one or more gene identifiers, separated by **tabs**. For example (shown with spacing):

```text
pathway_A    na    GENE1    GENE2    GENE3
regulon_B    na    GENE4    GENE5
```

Use actual tabs when creating this file. Descriptions are ignored. Blank lines
are skipped; duplicate set names and empty gene fields are rejected. Omit
`targets` when supplying `gene_sets_gmt`. All sets in the file are tested.

Run the same `Rscript scripts/run_nb.R config.txt` command (or container command).
Each set's response is the **sum of raw member counts per cell**. The existing
within-donor receptor test, covariates, donor random intercept, spline option and
fit checks apply to this response. Reference depth, mitochondrial counts and
number of detected genes exclude **every member of that set and every input
receptor**. No per-member detection filter is applied; detection thresholds refer
to cells with a positive summed count. Overlapping sets are allowed and fitted
separately. Duplicate members within a set are counted once. Missing identifiers,
empty sets, duplicate set names and sets containing input receptors are rejected;
remove receptors explicitly when defining the intended target set.

This asks whether receptor expression is associated with the set's **total RNA
count relative to reference depth**, conditional on the model adjustments. It is
not a test that every member changes, an enrichment test against other genes, or
a signed/weighted GRN activity test. Highly expressed members can dominate and
opposing member changes can cancel. A linear rate ratio describes the summed
response. Gene-level follow-up can help explain a set result, but requires its own
prespecified multiple-testing family.

The NB model is an assumption for the aggregate response: sums of correlated NB
counts are not generally NB distributed. Inspect model adequacy before interpreting
asymptotic p-values; these implementation tests do not establish calibration.
The NB variance models are described in the
[glmmTMB family documentation](https://glmmtmb.github.io/glmmTMB/reference/nbinom2.html).

Set runs write `gene_set_results.csv` (including set size) instead of
`gene_results.csv`, plus `gene_set_members.csv` with the exact analyzed membership.
`receptor_effects.csv` uses a `gene_set` identifier and `analysis.rds` records the
analysis type and membership. BH correction is across eligible **sets in that
run**, including eligible failed fits in its denominator. Overlap makes tests
dependent; the existing BH assumptions and cross-run multiplicity caveats apply.

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
- Multiple receptors via the legacy `receptors` option or unpooled R API: a joint omnibus test conditional on all their donor-mean
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
