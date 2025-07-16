#!/usr/bin/env bash
# Activate the conda environment and run the chemistry extractor in testing mode

cd "$(dirname "$0")"
# Ensure conda is available in this shell
if [ -z "$CONDA_EXE" ]; then
    if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
        source "$HOME/miniconda3/etc/profile.d/conda.sh"
    elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
        source "$HOME/anaconda3/etc/profile.d/conda.sh"
    else
        echo "Could not find conda.sh. Please ensure conda is installed and available."
        exit 1
    fi
fi

conda activate chem_env
cd src
python3 -m domains.chemistry.main --env testing 

conda deactivate