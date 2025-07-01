# Config Reconciliation Script Design and Usage

## Purpose and Design Principles
- **Validation and reconciliation**: Ensures all YAML config files for the chemistry domain are present, consistent, and valid before downstream processing.
- **Single entry point**: Centralizes config validation logic for maintainability and reproducibility.
- **Environment-driven**: Loads all config file paths from environment variables (via EnvironmentLoader), supporting flexible deployment and testing.
- **Fail-fast**: Exits immediately with a clear error if any required config file is missing.

## What Config Files Are Reconciled
- Ontology YAML
- Entity config YAML
- Validation config YAML
- Source mapping YAML
- Conflict resolution YAML
- Extraction config YAML
- Metric units YAML

## How Config Paths Are Loaded
- Uses `EnvironmentLoader` to load all config paths from environment variables.
- By default, prints only the YAML config file paths being reconciled.
- If run with `PRINT_ENV=1` or `--print-env`, prints all environment variables loaded.

## Usage Instructions
- **Default:**
  ```bash
  python3 common_dictionary/src/config_reconciliation.py
  ```
  Prints only the YAML config file paths being reconciled.
- **Print all environment variables:**
  ```bash
  PRINT_ENV=1 python3 common_dictionary/src/config_reconciliation.py
  # or
  python3 common_dictionary/src/config_reconciliation.py --print-env
  ```

## Output and Error Handling
- Prints a summary of the ontology, entity mapping, and attribute reconciliation.
- Fails fast with a clear error if any required config path is missing.
- Prints detailed validation and reconciliation results for all entities and attributes.

## Extensibility Notes
- To add new config files, add their paths to the environment file and update the script to check for them.
- The script can be extended to perform deeper validation, cross-file checks, or to output machine-readable reports.
- The script is portable and can be run from any directory as long as the environment is set up correctly. 