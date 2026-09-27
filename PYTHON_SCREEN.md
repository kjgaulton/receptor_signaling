# Receptor expression prediction screen

Does measured receptor RNA expression improve prediction of another gene's
expression, beyond measured technical and biological covariates, in one cell type?

This initial implementation provides an **exploratory genome-wide prediction
ranking**, not a calibrated hypothesis test. It does not produce association
p-values, FDR calls, causal effects, or evidence of receptor activation.

## Model and estimand

For each eligible target gene, compare two histogram gradient-boosted regressors:

- Baseline: log reference UMI depth, log number of detected reference genes,
  mitochondrial fraction, and user-specified metadata covariates.
- Full: baseline features plus normalized expression of the supplied receptor(s).

Targets and receptors use `log1p(10000 * count / reference_depth)`. Squared-error
prediction does not require a Gaussian response distribution, but it emphasizes
the conditional mean on this transformed scale. It is not a count likelihood
or a model of absolute RNA abundance. Trees can capture nonlinear relationships
and interactions. Zeros remain zeros; no imputation or smoothing is applied.

All cells from a donor are held out together using GroupKFold. Each donor has equal
total training weight. Evaluation averages squared errors within each donor and
then across donors, so donors with more captured cells do not dominate the score.
Donor IDs must represent biological individuals, including across repeated samples.
Donor-held-out performance can still reflect donor-level biology or condition;
it does not establish a within-donor association.

The ranking score is `1 - full_MSE / baseline_MSE`. A value of 0.10 means a 10%
reduction in held-out, donor-balanced MSE; a negative value means worse prediction.
This is not ordinary R-squared. Donor consistency is descriptive, not a p-value.

For each target, reference depth, detected-gene count, and mitochondrial fraction
exclude **all input receptors and that target**. This avoids their direct inclusion
in nuisance features and normalization denominators. Shared normalization and
compositional effects can still remain; this does not eliminate all technical bias.
Inputs must contain the whole measured transcriptome, not only highly variable genes.

## Run

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -e .
.venv/bin/receptor-screen data.h5ad \
  --receptors IL7R \
  --cell-type 'CD4 T' --cell-type-key cell_type \
  --donor-key donor_id --layer counts \
  --covariates batch condition \
  --out results/il7r_cd4
```

Supply multiple receptors after `--receptors` for a **joint** comparison. To measure
each receptor individually, run separately. Unique contributions conditional on
other receptors require leave-one-receptor-out comparisons, which are not yet
implemented. Joint improvement alone does not establish synergy.

- `.var_names` must contain unique matching gene identifiers. Ensembl IDs work if
  receptors are specified as Ensembl IDs and a mitochondrial mask is provided.
- Default counts layer is `counts`; `--layer X` explicitly selects `.X`.
- Default mitochondrial prefix is `MT-`; use `--mito-prefix mt-` for appropriate
  mouse annotations or `--mito-key is_mito` for a boolean `.var` column.
- Five folds require at least five donors. Use `--folds` to change this (minimum
  two). Very few donors provide little evidence of generalization.
- Additional numeric/categorical `.obs` covariates must be complete. Categories
  are encoded inside training folds; unseen categories become all-zero indicators.
  Do not add existing total-count/mitochondrial/detected-gene QC columns: the
  pipeline recalculates them. Do not add receptor, target, or outcome-derived
  annotations as covariates.
- Cells should already pass dataset-appropriate doublet and low-quality filtering.
  The tool validates counts but does not perform sample QC.
- By default, targets need detection in at least 20 cells and 1% of selected cells.
  Filtered genes remain in the output with a status. Receptors are not target genes.
  Target detection filtering uses all selected cells, so this is a descriptive
  screen, not an untouched confirmatory test set.

Outputs: `gene_scores.csv`, per-gene `donor_scores.csv`, cell-level `folds.csv`, and
`run.json` with parameters and dependency versions. Existing output directories
are not overwritten. The input matrix stays sparse, but the entire selected dataset
is held in memory. Fits are sequential across genes; `--threads` controls threads
per fit. Cost is roughly 2 × genes × folds boosted models. Pilot on a separate small
dataset first. There is no checkpoint/resume or distributed execution yet.

## Covariate decisions

Depth, mitochondrial fraction and detected-gene count are included by default.
Detected-gene count is correlated with depth and can also encode real biology;
mitochondrial fraction can similarly reflect state. Compare sensitivity analyses
before interpreting biology. This version fixes those three QC features; their
ablation is a planned extension.

Consider batch, condition, cell cycle, and finer cell states only when justified
by the question. Adjusting for a downstream state can remove the signal of interest.
Transcriptome-derived PCs or embeddings can contain the target/receptor signal;
do not casually use them as adjustment variables. Hold out a separate cohort for
stronger validation when possible.

## What is needed for association testing

The screen is a first stage. A defensible inferential extension should predefine
the null and validate its calibration before attaching p-values:

1. For a **predictive claim**, train/freeze models using discovery donors, then
   compare paired prediction losses on independent confirmation donors. Any
   donor-level inference must respect donor independence and adequate donor count.
   Do not run a naive t-test on overlapping cross-validation folds or on cells.
2. For a **conditional-association claim**, a conditional randomization procedure
   could model receptor expression given covariates, including detection and donor
   structure. Its validity depends on the conditional sampling model. Simple global
   or within-donor shuffling does not preserve receptor–QC relationships. Coarse
   QC-bin shuffling is an approximation, not an exact solution.
3. Validate null calibration using simulated depth-only associations, donor effects,
   dropout, cell-state mixtures, and correlated receptors. Assess power for nonlinear
   signals and sparse receptors. Calibrate before applying multiplicity correction.
4. Apply an appropriate multiple-testing procedure to valid p-values across the
   prespecified receptor × gene × cell-type family. BH needs suitable dependence
   assumptions; no FDR claim is made by this implementation.

A count-based alternative is boosted Poisson deviance with a genuine log-depth
offset, or a negative-binomial/GAM sensitivity analysis. Merely adding depth as a
feature is not equivalent to an offset. A two-part model can separately assess
detection and positive expression, with its own prespecified testing family.
These alternatives are not implemented in this prototype.

## Validation

Run `.venv/bin/python -m unittest discover -s tests -v`. Synthetic tests check
donor separation, recovery of a nonlinear receptor signal, sparse/dense agreement,
and rejection of invalid counts. They establish implementation behavior, not
inferential calibration or biological validity.

## Method references

- [scikit-learn histogram gradient boosting](https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.HistGradientBoostingRegressor.html)
- [GroupKFold](https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.GroupKFold.html)
- [Zimmerman et al., pseudoreplication in single-cell studies](https://www.nature.com/articles/s41467-021-21038-1)
- [Squair et al., false discoveries in single-cell differential expression](https://www.nature.com/articles/s41467-021-25960-2)

The latter papers motivate respecting biological replication; they do not validate
this particular predictive algorithm.
