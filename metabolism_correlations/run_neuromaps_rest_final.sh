#!/bin/bash
# Runs neuromaps_rest_final.py in the "neuromaps_env" conda environment.
# Creates the environment (and installs neuromaps/netneurotools) on first use.
set -euo pipefail

CONDA_BASE=/export02/data/vikramn/conda
ENV_NAME=neuromaps_env
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source "${CONDA_BASE}/etc/profile.d/conda.sh"

if ! conda env list | grep -qE "^${ENV_NAME}[[:space:]]"; then
    echo "Creating conda environment '${ENV_NAME}'..."
    conda create -n "${ENV_NAME}" python=3.10 -y
    conda activate "${ENV_NAME}"
    pip install pandas numpy matplotlib scipy neuromaps netneurotools
else
    conda activate "${ENV_NAME}"
fi

python "${SCRIPT_DIR}/neuromaps_rest_final.py" "$@"
