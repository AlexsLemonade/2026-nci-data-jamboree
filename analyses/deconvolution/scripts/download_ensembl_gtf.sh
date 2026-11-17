#!/bin/bash
# Download the Ensembl GRCh38 gene annotation (GTF) used to build the gene
# category table (see make_gene_categories.R). Ensembl 104 matches the
# annotation used by ScPCA.
#
# Writes data/genome-references/Homo_sapiens.GRCh38.<release>.gtf.gz at the
# repository root (git-ignored) and verifies it against Ensembl's published
# checksum.
#
# Usage: bash download_ensembl_gtf.sh [release]   (default release: 104)

set -euo pipefail

release="${1:-104}"
gtf="Homo_sapiens.GRCh38.${release}.gtf.gz"
base_url="https://ftp.ensembl.org/pub/release-${release}/gtf/homo_sapiens"

repo_root="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
out_dir="${repo_root}/data/genome-references"
mkdir -p "${out_dir}"

if [[ -f "${out_dir}/${gtf}" ]]; then
  echo "Already present: ${out_dir}/${gtf}"
else
  echo "Downloading ${base_url}/${gtf}"
  curl -fsSL -o "${out_dir}/${gtf}.tmp" "${base_url}/${gtf}"
  mv "${out_dir}/${gtf}.tmp" "${out_dir}/${gtf}"
fi

# Ensembl CHECKSUMS files use the BSD `sum` algorithm
expected="$(curl -fsSL "${base_url}/CHECKSUMS" | awk -v f="${gtf}" '$3 == f {print $1, $2}')"
observed="$(sum "${out_dir}/${gtf}" | awk '{print $1, $2}')"
if [[ "${expected}" != "${observed}" ]]; then
  echo "Checksum mismatch for ${gtf}: expected '${expected}', got '${observed}'" >&2
  exit 1
fi
echo "Checksum OK (${observed}): ${out_dir}/${gtf}"
