import argparse
from dataclasses import asdict
from importlib.metadata import version
import json
from pathlib import Path

import anndata as ad

from .core import Settings, screen


def main():
    p = argparse.ArgumentParser(description="Exploratory donor-held-out receptor prediction screen (no association p-values)")
    p.add_argument("input", type=Path)
    p.add_argument("--receptors", nargs="+", required=True)
    p.add_argument("--cell-type", required=True)
    p.add_argument("--cell-type-key", default="cell_type")
    p.add_argument("--donor-key", required=True)
    p.add_argument("--layer", default="counts", help="Raw counts layer; use X explicitly for adata.X")
    p.add_argument("--mito-prefix", default="MT-", help="Case-sensitive gene-ID prefix; mouse often mt-")
    p.add_argument("--mito-key", help="Boolean column in adata.var; overrides prefix")
    p.add_argument("--covariates", nargs="*", default=[])
    p.add_argument("--folds", type=int, default=5)
    p.add_argument("--min-detected", type=int, default=20)
    p.add_argument("--min-fraction", type=float, default=0.01)
    p.add_argument("--max-iter", type=int, default=100)
    p.add_argument("--min-samples-leaf", type=int, default=30)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    if args.out.exists():
        p.error("Output directory already exists; choose a new path")
    a = ad.read_h5ad(args.input)
    a = a[a.obs[args.cell_type_key].eq(args.cell_type)].copy()
    if a.n_obs == 0:
        p.error("No cells match the requested cell type")
    counts = a.X if args.layer == "X" else a.layers[args.layer]
    if args.mito_key:
        mask = a.var[args.mito_key]
        if mask.isna().any() or mask.dtype != bool:
            p.error("--mito-key must name a nonmissing boolean column")
        mito = mask.to_numpy()
    else:
        mito = a.var_names.str.startswith(args.mito_prefix)
    settings = Settings(folds=args.folds, min_detected=args.min_detected,
                        min_fraction=args.min_fraction, max_iter=args.max_iter,
                        min_samples_leaf=args.min_samples_leaf, threads=args.threads)
    scores, donors, folds = screen(counts, a.var_names, a.obs, args.receptors,
                                    args.donor_key, mito, args.covariates, settings)
    args.out.mkdir(parents=True)
    scores.to_csv(args.out / "gene_scores.csv", index=False)
    donors.to_csv(args.out / "donor_scores.csv", index=False)
    folds.to_csv(args.out / "folds.csv", index=False)
    metadata = dict(arguments={k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
                    settings=asdict(settings), n_cells=a.n_obs, n_genes=a.n_vars,
                    versions={x: version(x) for x in ["numpy", "pandas", "scipy", "scikit-learn", "anndata"]},
                    interpretation="Exploratory predictive ranking; no calibrated p-values or FDR")
    (args.out / "run.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Wrote {len(scores)} gene records to {args.out}")


if __name__ == "__main__":
    main()
