# test_config1: Chemistry Extractor POC Configs

This directory contains all configuration files for the Extractor Proof-of-Concept (POC) for the chemistry domain, as described in [ExtractorPOCImplementation.md](../../../documentation/ExtractorPOCImplementation.md).

## Scope
- **Entities:** Only `Compound`, `Polymer`, and `Biopolymer` (for hierarchy demonstration)
- **Sources:** Only PubChem and RDKit
- **Purpose:** Minimal, realistic, and testable configs for extractor and review pipeline POC
- **Alignment:** All files strictly follow the structure and required fields of the templates in `env_templates/chemistry/config_templates/`

## Included Config Files
- `ontology.yaml` — Entity hierarchy and relationships for POC
- `entity_config.yaml` — Attributes and enrichment rules for POC entities
- `extraction_config.yaml` — Extraction sources and parameters (PubChem, RDKit)
- `conflict_resolution.yaml` — Conflict resolution rules for POC attributes
- `relationships_config.yaml` — Relationships relevant to POC entities
- `validation_config.yaml` — Validation rules for POC attributes
- `source_mapping.yaml` — Attribute-to-source mapping for PubChem and RDKit

## Usage
These configs are designed for use with the extractor and review pipeline as described in ExtractorPOCImplementation.md. They are minimal, focused, and ready for immediate testing and extension. 