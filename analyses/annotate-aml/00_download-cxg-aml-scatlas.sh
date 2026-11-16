#!/bin/bash
# Download AML scAtlas from CZI CELLxGENE

set -euo pipefail

# Run from the directory that this shell script is saved in
cd $(dirname "$0")
module_dir=$(pwd)

# Create data directory if it doesn't exist
data_dir="../../data/AML-scAtlas"
mkdir -p "${data_dir}"

# Download the AML scAtlas dataset from CZI CELLxGENE
wget -O "${data_dir}/AML_scAtlas.h5ad" https://datasets.cellxgene.cziscience.com/be1f3ced-6b8f-47b0-a3e2-7edc8d466abd.h5ad
