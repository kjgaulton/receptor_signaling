"""Sparse-count input; only a target and the small feature matrix are densified."""
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import OneHotEncoder
from threadpoolctl import threadpool_limits


@dataclass
class Settings:
    folds: int = 5
    min_detected: int = 20
    min_fraction: float = 0.01
    max_iter: int = 100
    min_samples_leaf: int = 30
    seed: int = 42
    threads: int = 1


def screen(counts, genes, obs, receptors, donor_key, mito, covariates=(), settings=None):
    """Return gene scores, donor scores, and fold assignments.

    Receptors are tested jointly. Score = 1 - full MSE / baseline MSE,
    where MSE is averaged within donor then equally across donors.
    """
    s = settings or Settings()
    if (s.folds < 2 or s.min_detected < 1 or not 0 <= s.min_fraction <= 1
            or s.max_iter < 1 or s.min_samples_leaf < 1 or s.threads < 1):
        raise ValueError("Invalid settings")
    genes = pd.Index(genes)
    if not genes.is_unique:
        raise ValueError("Gene identifiers must be unique; resolve duplicates explicitly")
    receptors = list(dict.fromkeys(receptors))
    if not receptors or not set(receptors).issubset(genes):
        raise ValueError("Provide at least one receptor; every receptor must match a gene ID")
    c = sparse.csc_matrix(counts, dtype=float)
    c.sum_duplicates()
    c.eliminate_zeros()
    if c.shape != (len(obs), len(genes)):
        raise ValueError("Counts shape does not match metadata")
    if not np.isfinite(c.data).all() or (c.data < 0).any() or not np.allclose(c.data, np.round(c.data)):
        raise ValueError("Supply finite, nonnegative raw integer UMI counts")
    mito = np.asarray(mito, dtype=bool)
    if mito.shape != (len(genes),) or not mito.any():
        raise ValueError("Mitochondrial mask must identify at least one gene")
    if obs[donor_key].isna().any():
        raise ValueError("Donor IDs must not be missing")
    groups = obs[donor_key].astype(str).to_numpy()
    if len(np.unique(groups)) < s.folds:
        raise ValueError("Number of donors must be at least the number of folds")
    if donor_key in covariates:
        raise ValueError("Donor ID is a splitting variable, not a predictor")
    extra = obs[list(covariates)].copy().reset_index(drop=True)
    if extra.isna().any().any():
        raise ValueError("Extra covariates contain missing values; resolve before analysis")
    numeric = list(extra.select_dtypes(include=[np.number]).columns)
    if numeric and not np.isfinite(extra[numeric].to_numpy()).all():
        raise ValueError("Numeric covariates must be finite")
    categorical = [x for x in extra if x not in numeric]
    for col in categorical:
        extra[col] = extra[col].astype(str)
    folds = list(GroupKFold(s.folds).split(np.zeros(len(obs)), groups=groups))
    fold_ids = np.full(len(obs), -1)
    encoded = []
    for f, (train, test) in enumerate(folds):
        fold_ids[test] = f
        if extra.shape[1]:
            encoder = ColumnTransformer([
                ("numeric", "passthrough", numeric),
                ("categorical", OneHotEncoder(handle_unknown="ignore", sparse_output=False), categorical)
            ])
            encoded.append((encoder.fit_transform(extra.iloc[train]), encoder.transform(extra.iloc[test])))
        else:
            encoded.append((np.empty((len(train), 0)), np.empty((len(test), 0))))
    ri = genes.get_indexer(receptors)
    r = c[:, ri].toarray()
    for k, name in enumerate(receptors):
        if np.count_nonzero(r[:, k]) < s.min_detected or np.ptp(r[:, k]) == 0:
            raise ValueError(f"Receptor {name} has insufficient expression/variation")
    # Reference totals omit receptors for all targets; omit the target below.
    total = np.asarray(c.sum(axis=1)).ravel() - r.sum(axis=1)
    detected = np.asarray((c > 0).sum(axis=1)).ravel() - (r > 0).sum(axis=1)
    mt = np.asarray(c[:, mito].sum(axis=1)).ravel() - r[:, mito[ri]].sum(axis=1)
    rows, donor_rows = [], []
    for j, gene in enumerate(genes):
        if gene in receptors:
            continue
        ycount = c[:, j].toarray().ravel()
        n_detected = int(np.count_nonzero(ycount))
        if n_detected < max(s.min_detected, s.min_fraction * len(obs)):
            rows.append(dict(gene=gene, status="low_detection", n_detected=n_detected))
            continue
        depth = total - ycount
        if (depth <= 0).any():
            raise ValueError(f"Zero reference depth for target {gene}; filter empty cells or supply whole-transcriptome counts")
        qc = np.column_stack([np.log1p(depth), np.log1p(detected - (ycount > 0)),
                              (mt - ycount * mito[j]) / depth])
        y = np.log1p(1e4 * ycount / depth)
        receptor_features = np.log1p(1e4 * r / depth[:, None])
        predictions = np.empty((len(obs), 2))
        with threadpool_limits(limits=s.threads):
            for f, (train, test) in enumerate(folds):
                xtrain = np.column_stack([qc[train], encoded[f][0]])
                xtest = np.column_stack([qc[test], encoded[f][1]])
                _, inv, sizes = np.unique(groups[train], return_inverse=True, return_counts=True)
                weights = 1.0 / sizes[inv]
                weights *= len(train) / weights.sum()
                for model_index in range(2):
                    a, b = xtrain, xtest
                    if model_index:
                        a = np.column_stack([a, receptor_features[train]])
                        b = np.column_stack([b, receptor_features[test]])
                    model = HistGradientBoostingRegressor(
                        loss="squared_error", max_iter=s.max_iter, max_leaf_nodes=7,
                        min_samples_leaf=s.min_samples_leaf, l2_regularization=1,
                        early_stopping=False, random_state=s.seed)
                    model.fit(a, y[train], sample_weight=weights)
                    predictions[test, model_index] = model.predict(b)
        errors = (y[:, None] - predictions) ** 2
        by_donor = []
        for donor in np.unique(groups):
            mask = groups == donor
            base, full = errors[mask].mean(axis=0)
            by_donor.append((base, full))
            donor_rows.append(dict(gene=gene, donor=donor, fold=int(fold_ids[mask][0]),
                                   n_cells=int(mask.sum()), baseline_mse=base, full_mse=full,
                                   improvement=base-full))
        base, full = np.mean(by_donor, axis=0)
        rows.append(dict(gene=gene, status="scored", n_detected=n_detected,
                         detection_fraction=n_detected/len(obs), baseline_mse=base,
                         full_mse=full, relative_mse_reduction=1-full/base if base > 0 else np.nan,
                         donor_fraction_improved=float(np.mean(np.diff(by_donor, axis=1) < 0))))
    result = pd.DataFrame(rows)
    if "relative_mse_reduction" in result:
        result = result.sort_values("relative_mse_reduction", ascending=False, na_position="last")
    assignments = pd.DataFrame(dict(cell=obs.index.astype(str), donor=groups, fold=fold_ids))
    return result, pd.DataFrame(donor_rows), assignments
