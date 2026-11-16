#!/usr/bin/env python3
"""Train an scVI -> scANVI reference on the AML scAtlas subset (donors aged 29 and under).

Steps (see .workspace/annotate-aml/reference.md, step 3):
  - restrict genes to those present in every library of the ScPCA AML project (SCPCP000007)
  - select highly variable genes on raw counts (`seurat_v3`, batch-aware by study)
  - train SCVI with `Sample` as the batch, using scArches-friendly architecture settings so
    ScPCA queries can be mapped later
  - train SCANVI from the SCVI model on `label_fine`
  - save both models, the HVG list, training histories, and the reference AnnData with its
    latent spaces

Query mapping later (scArches):
  SCANVI.prepare_query_anndata(query, model_dir); SCANVI.load_query_data(query, model_dir);
  then a short fine-tune with `train(max_epochs=..., plan_kwargs={"weight_decay": 0.0})`.

Run from anywhere with:
  uv run --project analyses/annotate-aml python analyses/annotate-aml/04_train-scanvi-ref.py

For validation (step 4), leave one study out; outputs go to `scanvi_ref_no-<study>/` so the
full reference isn't overwritten:
  uv run --project analyses/annotate-aml python analyses/annotate-aml/04_train-scanvi-ref.py \
    --holdout-study mumme_2023
"""

import argparse
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
from anndata.io import read_elem

# ---------------------------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------------------------

SCPCA_PROJECT = "SCPCP000007"  # ScPCA project made up primarily of AML samples
LABEL_KEY = "label_fine"
BATCH_KEY = "Sample"
HVG_BATCH_KEY = "Study"  # per-study HVG; many samples are too small for per-sample loess
UNLABELED_CATEGORY = "Unknown"  # required by SCANVI; every reference cell is labeled

# scArches-friendly architecture (scvi-tools reference mapping tutorial)
SCVI_KWARGS = dict(
    n_latent=30,
    n_layers=2,
    gene_likelihood="nb",
    dropout_rate=0.2,
    use_layer_norm="both",
    use_batch_norm="none",
    encode_covariates=True,
)


def find_root(start: Path) -> Path:
    """Return the first parent directory of `start` that contains a .git directory."""
    for path in [start, *start.parents]:
        if (path / ".git").exists():
            return path
    raise FileNotFoundError(f"No .git directory found above {start}")


def scpca_shared_genes(project_dir: Path) -> pd.Index:
    """Ensembl IDs present in every `*_filtered_rna.h5ad` library under `project_dir`."""
    files = sorted(project_dir.rglob("*_filtered_rna.h5ad"))
    if not files:
        raise FileNotFoundError(f"No *_filtered_rna.h5ad files under {project_dir}")
    shared = None
    for path in files:
        with h5py.File(path, "r") as f:
            genes = pd.Index(read_elem(f["var"]).index)
        shared = genes if shared is None else shared.intersection(genes)
    print(f"{len(files)} ScPCA libraries in {project_dir.name}; {len(shared):,} shared genes")
    return shared


def main() -> None:
    root = find_root(Path(__file__).resolve().parent)
    module_dir = root / "analyses" / "annotate-aml"

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--input", type=Path, default=root / "data" / "AML-scAtlas" / "AML_scAtlas_le29.h5ad"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="default: results/scanvi_ref (or results/scanvi_ref_no-<study> with a holdout)",
    )
    parser.add_argument(
        "--holdout-study",
        default=None,
        help="leave this Study out of training (removed before gene selection)",
    )
    parser.add_argument("--n-hvg", type=int, default=3000)
    parser.add_argument(
        "--scvi-max-epochs",
        type=int,
        default=None,
        help="default: scvi-tools heuristic based on the number of cells (with early stopping)",
    )
    parser.add_argument("--scanvi-max-epochs", type=int, default=20)
    parser.add_argument("--n-samples-per-label", type=int, default=100)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args()

    if args.output_dir is None:
        name = "scanvi_ref" + (f"_no-{args.holdout_study}" if args.holdout_study else "")
        args.output_dir = module_dir / "results" / name

    scvi.settings.seed = args.seed
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # --- data ---------------------------------------------------------------------------------
    adata = ad.read_h5ad(args.input)
    print(f"Input: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    if not np.array_equal(adata.X.data, np.round(adata.X.data)):
        raise ValueError("`X` must hold raw integer counts")

    if args.holdout_study:
        is_holdout = (adata.obs["Study"] == args.holdout_study).to_numpy()
        if not is_holdout.any():
            raise ValueError(f"Study {args.holdout_study!r} not found")
        print(f"Holding out {is_holdout.sum():,} cells from {args.holdout_study}")
        adata = adata[~is_holdout].copy()
        for col in adata.obs.columns:
            if isinstance(adata.obs[col].dtype, pd.CategoricalDtype):
                adata.obs[col] = adata.obs[col].cat.remove_unused_categories()
        print(f"Training on {adata.n_obs:,} cells")

    # --- genes --------------------------------------------------------------------------------
    shared = scpca_shared_genes(root / "data" / "scpca" / SCPCA_PROJECT)
    keep_genes = adata.var_names.isin(shared)
    print(f"Atlas genes shared with ScPCA: {keep_genes.sum():,} of {adata.n_vars:,}")
    if keep_genes.sum() < 0.9 * adata.n_vars:
        raise ValueError("Fewer than 90% of atlas genes are shared with ScPCA")
    adata = adata[:, keep_genes].copy()

    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=args.n_hvg,
        batch_key=HVG_BATCH_KEY,
        subset=True,
    )
    print(f"Selected {adata.n_vars:,} HVGs (seurat_v3, batch_key={HVG_BATCH_KEY})")
    adata.var.to_csv(args.output_dir / "hvg.tsv", sep="\t")

    # --- SCVI ---------------------------------------------------------------------------------
    scvi.model.SCVI.setup_anndata(adata, batch_key=BATCH_KEY)
    scvi_model = scvi.model.SCVI(adata, **SCVI_KWARGS)
    print(scvi_model)
    scvi_model.train(max_epochs=args.scvi_max_epochs, early_stopping=True)
    scvi_model.save(args.output_dir / "scvi_model", overwrite=True, save_anndata=False)
    pd.concat(scvi_model.history, axis=1).to_csv(args.output_dir / "scvi_history.tsv", sep="\t")
    adata.obsm["X_scVI"] = scvi_model.get_latent_representation()

    # --- SCANVI -------------------------------------------------------------------------------
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        unlabeled_category=UNLABELED_CATEGORY,
        labels_key=LABEL_KEY,
    )
    scanvi_model.train(
        max_epochs=args.scanvi_max_epochs,
        n_samples_per_label=args.n_samples_per_label,
    )
    scanvi_model.save(args.output_dir / "scanvi_model", overwrite=True, save_anndata=False)
    pd.concat(scanvi_model.history, axis=1).to_csv(
        args.output_dir / "scanvi_history.tsv", sep="\t"
    )
    adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation()

    # --- sanity check: self-prediction on the training data (optimistic; real validation is
    # step 4) -------------------------------------------------------------------------------
    adata.obs["scanvi_pred"] = pd.Categorical(scanvi_model.predict())
    correct = adata.obs["scanvi_pred"].astype(str) == adata.obs[LABEL_KEY].astype(str)
    per_label = (
        adata.obs.assign(correct=correct)
        .groupby(LABEL_KEY, observed=True)["correct"]
        .agg(["size", "mean"])
        .rename(columns={"size": "n_cells", "mean": "self_accuracy"})
        .sort_values("self_accuracy")
    )
    print(f"\nOverall self-prediction accuracy: {correct.mean():.3f}")
    print(per_label.to_string(float_format="{:.3f}".format))
    per_label.to_csv(args.output_dir / "scanvi_self_accuracy.tsv", sep="\t")

    # --- save reference AnnData ---------------------------------------------------------------
    adata.uns["scanvi_ref_info"] = {
        "scpca_project": SCPCA_PROJECT,
        "n_hvg": args.n_hvg,
        "hvg_batch_key": HVG_BATCH_KEY,
        "batch_key": BATCH_KEY,
        "labels_key": LABEL_KEY,
        "seed": args.seed,
        "holdout_study": args.holdout_study or "",
    }
    ref_file = args.output_dir / "AML_scAtlas_le29_scanvi_ref.h5ad"
    adata.write_h5ad(ref_file, compression="gzip")
    print(f"\nWrote models and reference to {args.output_dir}")


if __name__ == "__main__":
    main()
