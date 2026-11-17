#!/usr/bin/env python3
"""Subset the AML scAtlas to donors aged 29 and under and slim it down for reference building.

Steps (see .workspace/annotate-aml/reference.md, step 1):
  - keep cells with numeric `age_floor <= 29` (unknown ages are excluded)
  - drop Smart-seq2 cells
  - make the raw counts (`raw/X`) the main matrix
  - keep a reduced set of metadata columns
  - add `label_broad`, `label_fine`, and `malignancy` columns
  - drop cells in `label_fine` classes with fewer than `--min-cells` cells

The input file is never loaded into memory in full: obs is read on its own, and only the selected
rows of `raw/X` are read from disk.

Run from anywhere with:
  uv run --project analyses/annotate-aml python analyses/annotate-aml/01_subset-atlas.py
"""

import argparse
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, sparse_dataset

# ---------------------------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------------------------

MAX_AGE = 29
EXCLUDED_ASSAYS = ["Smart-seq2"]
HEALTHY_DISEASE = "normal"
HSPC_CELL_TYPE = "hematopoietic stem cell"

# obs columns carried into the output
KEEP_OBS = [
    "donor_id",
    "Sample",
    "Study",
    "assay",
    "assay_ontology_term_id",
    "disease",
    "disease_ontology_term_id",
    "age_floor",
    "development_stage",
    "ELN Classification",
    "ELN Risk Group",
    "Cytogenetics",
    "Translocations",
    "is_primary_data",
    "cell_type",
    "cell_type_ontology_term_id",
    "Author Cell Type",
    "HSPC Cell Type",
]


def find_root(start: Path) -> Path:
    """Return the first parent directory of `start` that contains a .git directory."""
    for path in [start, *start.parents]:
        if (path / ".git").exists():
            return path
    raise FileNotFoundError(f"No .git directory found above {start}")


def build_labels(obs: pd.DataFrame) -> pd.DataFrame:
    """Add label_broad, label_fine, and malignancy columns to obs (returns a copy)."""
    obs = obs.copy()
    cell_type = obs["cell_type"].astype(str)
    # the source stores missing sub-labels as the literal string "NA" (R export), so treat
    # those as missing too
    hspc_sub = obs["HSPC Cell Type"].astype("string")
    hspc_sub = hspc_sub.mask(hspc_sub.isin(["NA", "nan", ""]))

    # label_broad: cell_type as is
    obs["label_broad"] = cell_type

    # label_fine: HSPC class replaced by its prefixed sub-label; HSPC cells with no sub-label
    # become "HSPC"
    is_hspc = cell_type == HSPC_CELL_TYPE
    fine = cell_type.copy()
    fine[is_hspc] = ("HSPC:" + hspc_sub[is_hspc]).fillna("HSPC").astype(str)
    obs["label_fine"] = fine

    # malignancy: `-like` HSPC sub-labels and LSC are malignant-like; every cell from a healthy
    # donor is normal (this takes precedence: a few `-like` cells come from healthy donors);
    # everything else is unknown
    is_malignant_like = (
        is_hspc & (hspc_sub.str.endswith("-like") | (hspc_sub == "LSC")).fillna(False)
    )
    is_healthy = obs["disease"].astype(str) == HEALTHY_DISEASE
    malignancy = pd.Series("unknown", index=obs.index)
    malignancy[is_malignant_like] = "malignant-like"
    malignancy[is_healthy] = "normal"
    n_overridden = int((is_malignant_like & is_healthy).sum())
    print(f"  {n_overridden} `-like`/LSC cells from healthy donors set to `normal`")
    obs["malignancy"] = malignancy

    for col in ["label_broad", "label_fine", "malignancy"]:
        obs[col] = obs[col].astype("category")
    return obs


def main() -> None:
    root = find_root(Path(__file__).resolve().parent)
    data_dir = root / "data" / "AML-scAtlas"

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", type=Path, default=data_dir / "AML_scAtlas.h5ad")
    parser.add_argument("--output", type=Path, default=data_dir / "AML_scAtlas_le29.h5ad")
    parser.add_argument(
        "--min-cells",
        type=int,
        default=100,
        help="drop label_fine classes with fewer than this many cells",
    )
    args = parser.parse_args()

    with h5py.File(args.input, "r") as f:
        # --- metadata -------------------------------------------------------------------------
        obs = read_elem(f["obs"])
        var = read_elem(f["var"])
        raw_var = read_elem(f["raw/var"])
        print(f"Input: {obs.shape[0]:,} cells x {var.shape[0]:,} genes")

        if not var.index.equals(raw_var.index):
            raise ValueError("`var` and `raw/var` have different gene indexes")

        missing = [c for c in KEEP_OBS if c not in obs.columns]
        if missing:
            raise KeyError(f"Expected obs columns not found: {missing}")

        # --- cell filters ---------------------------------------------------------------------
        age = pd.to_numeric(obs["age_floor"].astype(str), errors="coerce")
        keep_age = (age <= MAX_AGE).to_numpy()
        print(f"Age <= {MAX_AGE} (unknown excluded): {keep_age.sum():,} cells")

        is_excluded_assay = obs["assay"].astype(str).isin(EXCLUDED_ASSAYS).to_numpy()
        print(
            f"Dropping {(keep_age & is_excluded_assay).sum():,} cells from "
            f"{EXCLUDED_ASSAYS} assays"
        )
        keep = keep_age & ~is_excluded_assay

        # --- labels and tiny-class filter -----------------------------------------------------
        print("Building labels")
        sub_obs = build_labels(obs.loc[keep, KEEP_OBS])

        fine_counts = sub_obs["label_fine"].value_counts()
        small = fine_counts[fine_counts < args.min_cells]
        if len(small) > 0:
            dropped = sub_obs[sub_obs["label_fine"].isin(small.index)]
            print(
                f"Dropping {len(small)} label_fine classes with < {args.min_cells} cells "
                f"({small.sum():,} cells):"
            )
            print(
                pd.crosstab(dropped["label_fine"].astype(str), dropped["malignancy"].astype(str))
                .assign(total=lambda df: df.sum(axis=1))
                .sort_values("total", ascending=False)
                .to_string()
            )
        keep[keep] = ~sub_obs["label_fine"].isin(small.index).to_numpy()
        sub_obs = sub_obs[~sub_obs["label_fine"].isin(small.index)]

        # --- counts ---------------------------------------------------------------------------
        idx = np.flatnonzero(keep)
        print(f"Reading raw counts for {len(idx):,} cells")
        counts = sparse_dataset(f["raw/X"])[idx]

    # sanity check that the matrix really holds integer counts
    if not np.array_equal(counts.data, np.round(counts.data)):
        raise ValueError("`raw/X` contains non-integer values")

    # tidy categoricals so dropped levels don't linger
    for col in sub_obs.columns:
        if isinstance(sub_obs[col].dtype, pd.CategoricalDtype):
            sub_obs[col] = sub_obs[col].cat.remove_unused_categories()

    adata = ad.AnnData(X=counts, obs=sub_obs, var=var)
    adata.uns["subset_info"] = {
        "source": args.input.name,
        "max_age_floor": MAX_AGE,
        "excluded_assays": EXCLUDED_ASSAYS,
        "min_cells_label_fine": args.min_cells,
        "dropped_label_fine": list(small.index.astype(str)),
        "X": "raw counts from raw/X of the source file",
    }

    # --- summary ------------------------------------------------------------------------------
    print(f"\nOutput: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print(f"Donors: {adata.obs['donor_id'].nunique()}, Studies: {adata.obs['Study'].nunique()}")
    print("\nCells by disease:")
    print(adata.obs["disease"].value_counts().to_string())
    print("\nCells by assay:")
    print(adata.obs["assay"].value_counts().to_string())
    print("\nCells by malignancy:")
    print(adata.obs["malignancy"].value_counts().to_string())
    print("\nCells by label_fine:")
    print(adata.obs["label_fine"].value_counts().to_string())

    args.output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output, compression="gzip")
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
