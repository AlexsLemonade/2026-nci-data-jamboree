"""Split the AML scAtlas CELLxGENE h5ad into one h5ad per reference sample.

build_reference.R reads one input file at a time, so per-sample files keep its
peak memory to the largest sample instead of the whole atlas.

Each output file has:
  - X: raw integer counts (from the source's raw/X), cells x genes, all genes
    in the source order
  - var: Ensembl IDs as the index, plus a gene_symbol column
  - obs: donor_id (the reference sample), deconv_cell_type, deconv_cell_state
    and a few source columns for provenance

Labels:
  - deconv_cell_type is "Malignant" for cells from AML donors whose
    "HSPC Cell Type" ends in "-like" or is "LSC"; otherwise it is the
    "Author Cell Type".
  - deconv_cell_state is the "HSPC Cell Type" label for malignant cells, and
    the cell type for all others.

Cells from Smart-seq2 libraries are dropped: their counts are reads, not UMIs.

A manifest (samples.tsv) lists each written file with its cell counts.

Usage:
  uv run python scripts/split_by_sample.py [--input <h5ad>] [--output-dir <dir>]
"""

import argparse
import subprocess
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

MALIGNANT = "Malignant"
EXCLUDED_ASSAYS = {"Smart-seq2"}
AML_DISEASE = "myeloid leukemia"

# source obs column -> output obs column
KEEP_OBS = {
    "donor_id": "donor_id",
    "Sample": "library",
    "Study": "study",
    "assay": "assay",
    "disease": "disease",
    "Main Age Group": "main_age_group",
    "Author Cell Type": "author_cell_type",
    "HSPC Cell Type": "hspc_cell_type",
}


def repo_root() -> Path:
    out = subprocess.run(
        ["git", "-C", str(Path(__file__).parent), "rev-parse", "--show-toplevel"],
        check=True,
        capture_output=True,
        text=True,
    )
    return Path(out.stdout.strip())


def parse_args() -> argparse.Namespace:
    data_dir = repo_root() / "data" / "aml-scatlas"
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--input",
        type=Path,
        default=data_dir / "aml-scatlas_cellxgene.h5ad",
        help="CELLxGENE AML scAtlas h5ad",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=data_dir / "samples",
        help="directory for per-sample h5ad files and samples.tsv",
    )
    return parser.parse_args()


def assign_labels(obs: pd.DataFrame) -> pd.DataFrame:
    """Add deconv_cell_type and deconv_cell_state columns."""
    author = obs["author_cell_type"].astype(str)
    hspc = obs["hspc_cell_type"].astype("string").fillna("")
    is_malignant_label = hspc.str.endswith("-like") | (hspc == "LSC")
    malignant = ((obs["disease"] == AML_DISEASE) & is_malignant_label).to_numpy()
    obs["deconv_cell_type"] = np.where(malignant, MALIGNANT, author)
    obs["deconv_cell_state"] = np.where(malignant, hspc, author)
    return obs


def read_rows(x: h5py.Group, rows: np.ndarray, n_genes: int) -> sp.csr_matrix:
    """Read the given (sorted) rows of a CSR matrix stored in an h5ad group."""
    indptr = x["indptr"]
    # read each run of consecutive rows as one slice
    breaks = np.flatnonzero(np.diff(rows) != 1) + 1
    data, indices, row_nnz = [], [], []
    for run in np.split(rows, breaks):
        start, end = int(indptr[run[0]]), int(indptr[run[-1] + 1])
        run_ptr = indptr[run[0] : run[-1] + 2]
        data.append(x["data"][start:end])
        indices.append(x["indices"][start:end])
        row_nnz.append(np.diff(run_ptr))
    data = np.concatenate(data)
    if not np.array_equal(data, np.round(data)):
        raise ValueError("raw/X contains non-integer values")
    new_ptr = np.concatenate([[0], np.cumsum(np.concatenate(row_nnz))])
    return sp.csr_matrix(
        (data.astype(np.int32), np.concatenate(indices), new_ptr),
        shape=(len(rows), n_genes),
    )


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # cells of one donor are scattered through the file, so cache many chunks
    with h5py.File(args.input, "r", rdcc_nbytes=512 * 2**20, rdcc_nslots=100_003) as f:
        obs = ad.io.read_elem(f["obs"])
        raw_var = ad.io.read_elem(f["raw/var"])
        n_genes = raw_var.shape[0]

        obs = obs[list(KEEP_OBS)].rename(columns=KEEP_OBS)
        obs = assign_labels(obs)
        var = pd.DataFrame(
            {"gene_symbol": raw_var["gene_symbol"].astype(str).to_numpy()},
            index=pd.Index(raw_var.index.astype(str), name=None),
        )

        keep = ~obs["assay"].isin(EXCLUDED_ASSAYS).to_numpy()
        print(f"Dropping {int((~keep).sum())} cells from {sorted(EXCLUDED_ASSAYS)}")

        manifest = []
        donors = obs.loc[keep, "donor_id"].astype(str)
        for donor in sorted(donors.unique()):
            rows = np.flatnonzero(
                keep & (obs["donor_id"].astype(str) == donor).to_numpy()
            )
            sample_obs = obs.iloc[rows].copy()
            for column in sample_obs.select_dtypes("category").columns:
                sample_obs[column] = sample_obs[column].cat.remove_unused_categories()
            counts = read_rows(f["raw/X"], rows, n_genes)
            out = args.output_dir / f"{donor}.h5ad"
            ad.AnnData(X=counts, obs=sample_obs, var=var).write_h5ad(out)
            manifest.append(
                {
                    "donor_id": donor,
                    "file": out.name,
                    "n_cells": len(rows),
                    "n_malignant": int(
                        (sample_obs["deconv_cell_type"] == MALIGNANT).sum()
                    ),
                    "study": ";".join(sorted(sample_obs["study"].unique())),
                    "disease": ";".join(sorted(sample_obs["disease"].unique())),
                }
            )
            print(f"Wrote {out.name}: {len(rows)} cells")

    pd.DataFrame(manifest).to_csv(
        args.output_dir / "samples.tsv", sep="\t", index=False
    )
    print(f"Wrote {len(manifest)} samples to {args.output_dir}")


if __name__ == "__main__":
    main()
