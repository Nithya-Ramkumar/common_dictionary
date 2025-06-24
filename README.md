# Common Dictionary Module

A modular, extensible system for extracting, reconciling, validating, enriching, and reviewing entities and relationships in scientific domains (starting with chemistry). Emphasizes provenance, traceability, human-in-the-loop review, and standards-compliant outputs for downstream applications such as knowledge graph construction and embedding-based search.

## Setup

1. Clone the repository and navigate to the `common_dictionary` directory.
2. Create a Python virtual environment:
   ```bash
   python3 -m venv venv
   source venv/bin/activate
   ```
3. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

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