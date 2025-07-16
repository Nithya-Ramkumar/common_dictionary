# Common Dictionary Module

A modular, extensible system for extracting, reconciling, validating, enriching, and reviewing entities and relationships in scientific domains (starting with chemistry). Emphasizes provenance, traceability, human-in-the-loop review, and standards-compliant outputs for downstream applications such as knowledge graph construction and embedding-based search.

## Setup (Recommended: Conda Environment)

1. Clone the repository and navigate to the `common_dictionary` directory.
2. Create the conda environment (recommended for RDKit and scientific packages):
   ```bash
   conda env create -f environment.yml
   conda activate chem_env
   ```
   If you have not yet created `environment.yml`, you can export your current environment with:
   ```bash
   conda env export -n chem_env > environment.yml
   ```
3. Install any additional dependencies (if needed):
   ```bash
   pip install -r requirements.txt
   ```
4. (Optional) Remove the old Python venv to free up space:
   ```bash
   rm -rf .venv
   ```

## Migration Note
- The project has migrated from a Python venv-based setup to a conda environment for better compatibility with RDKit and scientific libraries.
- All environment management should now be done with conda.
- The `.venv` directory is no longer needed and can be safely deleted.

## Modules

### config_reconciliation.py
Validates that all entity attributes in `entity_config.yaml` are mapped in `source_priority.yaml` and `conflict_resolution.yaml`, and that there are no orphaned mappings.

### extractor.py
*Coming soon...*

### reconciler.py
*Coming soon...*

### ...

## Usage

Instructions for running each module will be added as they are implemented. 