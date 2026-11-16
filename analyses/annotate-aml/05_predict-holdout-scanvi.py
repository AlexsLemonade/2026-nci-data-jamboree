#!/usr/bin/env python3
"""Map held-out atlas cells onto the scANVI reference trained without them (scArches) and predict labels.

See .workspace/annotate-aml/evaluation.md (design step 2). The held-out cells are treated
exactly like a new query (e.g., ScPCA):
  - their labels are hidden (set to the unlabeled category), so nothing leaks from the truth
  - each `Sample` is a new batch
  - SCANVI.prepare_query_anndata -> SCANVI.load_query_data -> short fine-tune -> predict

Input:
  - data/AML-scAtlas/AML_scAtlas_le29.h5ad (02_subset-atlas.py)
  - results/scanvi_ref_no-<study>/scanvi_model/ (04_train-scanvi-ref.py --holdout-study <study>)
Output (results/holdout/):
  - scanvi_predictions_<study>.tsv: cell_id, scanvi_label, scanvi_max_prob
  - scanvi_probabilities_<study>.tsv.gz: per-label probabilities
  - scanvi_query_<study>.h5ad: obs (cell_id, Sample) + obsm["X_scANVI"], no expression matrix
  - scanvi_query_model_<study>/: the fine-tuned query model

Run from anywhere with:
  uv run --project analyses/annotate-aml python analyses/annotate-aml/05_predict-holdout-scanvi.py
"""

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scvi

LABEL_KEY = "label_fine"
BATCH_KEY = "Sample"


def find_root(start: Path) -> Path:
    """Return the first parent directory of `start` that contains a .git directory."""
    for path in [start, *start.parents]:
        if (path / ".git").exists():
            return path
    raise FileNotFoundError(f"No .git directory found above {start}")


def main() -> None:
    root = find_root(Path(__file__).resolve().parent)
    module_dir = root / "analyses" / "annotate-aml"

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--holdout-study", default="mumme_2023")
    parser.add_argument(
        "--input", type=Path, default=root / "data" / "AML-scAtlas" / "AML_scAtlas_le29.h5ad"
    )
    parser.add_argument(
        "--model-dir",
        type=Path,
        default=None,
        help="default: results/scanvi_ref_no-<study>/scanvi_model",
    )
    parser.add_argument("--output-dir", type=Path, default=module_dir / "results" / "holdout")
    parser.add_argument("--max-epochs", type=int, default=100)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args()

    study = args.holdout_study
    if args.model_dir is None:
        args.model_dir = module_dir / "results" / f"scanvi_ref_no-{study}" / "scanvi_model"
    if not args.model_dir.exists():
        raise FileNotFoundError(
            f"{args.model_dir} not found; run 04_train-scanvi-ref.py --holdout-study {study}"
        )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    scvi.settings.seed = args.seed

    # --- query cells --------------------------------------------------------------------------
    backed = ad.read_h5ad(args.input, backed="r")
    is_holdout = (backed.obs["Study"] == study).to_numpy()
    if not is_holdout.any():
        raise ValueError(f"Study {study!r} not found in {args.input}")
    query = backed[is_holdout].to_memory()
    backed.file.close()
    print(f"Query: {query.n_obs:,} cells from {study}")

    # keep only what the model needs, and hide the true labels. The model directory (not a
    # loaded model) is passed to scArches: loading it on its own would need the reference AnnData.
    reference_model = str(args.model_dir)
    registry = scvi.model.SCANVI.load_registry(reference_model)
    unlabeled = registry["field_registries"]["labels"]["state_registry"]["unlabeled_category"]
    query.obs = pd.DataFrame(
        {
            BATCH_KEY: pd.Categorical(query.obs[BATCH_KEY].astype(str)),
            LABEL_KEY: pd.Categorical([unlabeled] * query.n_obs),
        },
        index=query.obs_names,
    )
    print(f"{query.obs[BATCH_KEY].nunique()} query samples; labels set to {unlabeled!r}")

    # --- scArches -----------------------------------------------------------------------------
    query = scvi.model.SCANVI.prepare_query_anndata(query, reference_model, inplace=False)
    print(f"Query genes after matching the reference: {query.n_vars:,}")

    query_model = scvi.model.SCANVI.load_query_data(query, reference_model)
    query_model.train(
        max_epochs=args.max_epochs,
        plan_kwargs={"weight_decay": 0.0},
        check_val_every_n_epoch=10,
    )
    query_model.save(
        args.output_dir / f"scanvi_query_model_{study}", overwrite=True, save_anndata=False
    )

    # --- predictions --------------------------------------------------------------------------
    probs = query_model.predict(soft=True)
    probs.index.name = "cell_id"
    pred = pd.DataFrame(
        {
            "scanvi_label": probs.idxmax(axis=1),
            "scanvi_max_prob": probs.max(axis=1),
        },
        index=probs.index,
    )
    pred.to_csv(args.output_dir / f"scanvi_predictions_{study}.tsv", sep="\t")
    probs.to_csv(args.output_dir / f"scanvi_probabilities_{study}.tsv.gz", sep="\t")
    print("\nPredicted label counts:")
    print(pred["scanvi_label"].value_counts().to_string())
    print(f"\nMedian max probability: {np.median(pred['scanvi_max_prob']):.3f}")

    # --- latent space for plotting ------------------------------------------------------------
    latent = ad.AnnData(
        obs=query.obs[[BATCH_KEY]].copy(),
        obsm={"X_scANVI": query_model.get_latent_representation()},
    )
    latent.write_h5ad(args.output_dir / f"scanvi_query_{study}.h5ad")
    print(f"\nWrote predictions to {args.output_dir}")


if __name__ == "__main__":
    main()
