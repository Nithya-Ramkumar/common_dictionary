# Extractor Proof-of-Concept (POC) Implementation Plan

---

## 1. Overview

This document outlines the design and implementation plan for a focused Extractor POC in the Common Dictionary Module. The goal is to demonstrate a minimal, robust pipeline for extracting and enriching `Compound` and `Polymer` entities using only PubChem and RDKit, with strong provenance, UUID, and confidence tracking, and a user-friendly review interface.

---

## 2. Scope & Requirements

- **Entities:** Only `Compound` and `Polymer`
- **Sources:** Only PubChem and RDKit
- **Enrichment:** UUID and confidence score required in output
- **Config:** Minimal, test-focused config files in a dedicated directory
- **Environment:** All config/data paths set via an absolute-path environment file, loaded by a single environment loader
- **Output:** Schema-validated JSON for entities and relationships
- **UI:** Streamlit-based, easy-to-navigate review and edit interface

---

## 3. Directory & Config Structure

### 3.1. Config Directory
- `common_dictionary/config/test_config1/`
  - `ontology.yaml`: Defines only `Compound` and `Polymer` entities, their attributes, and allowed relationships
  - `entity_config.yaml`: Extraction/enrichment rules for each entity type
  - `extraction_config.yaml`: Only PubChem and RDKit as sources, with endpoints and logic

### 3.2. Output Schema Template
- `common_dictionary/env_templates/output_schemas.yaml`: JSON schema for entities and relationships, including UUID, confidence, provenance, etc.

### 3.3. Environment File
- `common_dictionary/env_templates/env.testing`: Contains **absolute paths** for all config and output directories, e.g.:
  ```env
  CONFIG_DIR=/home/nithya-ramkumar/Programming/LLM-Test-New/common_dictionary/config/test_config1
  DATA_OUTPUT_DIR=/home/nithya-ramkumar/Programming/LLM-Test-New/common_dictionary/data/outputs/testing
  OUTPUT_SCHEMA=/home/nithya-ramkumar/Programming/LLM-Test-New/common_dictionary/env_templates/output_schemas.yaml
  ```

---

## 4. Environment Loader

- Reads the `COMMON_DICT_ENV` environment variable (e.g., `testing`)
- Loads the corresponding `env.<env>` file from `env_templates/`
- All config and data paths are derived from this loader
- **All modules** (including `config_reconciliation.py` and `extractor.py`) use only this loader for config access

---

## 5. Extractor Architecture

### 5.1. Modular Source Extractors
- Each source (PubChem, RDKit) has its own extractor module/class
- Handles extraction/enrichment for its source, returns standardized dicts

### 5.2. Orchestrator (Main Extractor)
- Loads config via `env_loader`
- For each entity, calls the appropriate source extractors as per config
- Merges results, attaches provenance, generates UUIDs, computes confidence scores
- Validates output against `output_schemas.yaml`
- Writes `entities.json` and `relationships.json` (optionally prefixed with `raw_`) to the output directory

### 5.3. Extensibility
- New sources can be added as new extractor modules/classes
- Fallback logic can be added in the orchestrator layer in the future

---

## 6. Output Files

- `entities.json`: List of entities with all attributes, UUID, provenance, confidence, and review history (if any)
- `relationships.json`: List of relationships with UUID, from/to entity references, provenance, confidence, and review history
- All outputs go to the directory specified in the environment file

---

## 7. Streamlit UI

### 7.1. Navigation
- Sidebar: Switch between entities, relationships, and review logs
- Main area: Table or card view of entities/relationships, with search/filter
- Detail/edit view: Click to expand/edit an entity or relationship
- Review actions: Approve, reject, or edit attributes; all actions logged

### 7.2. Features
- Schema validation on edit
- Download/upload for outputs
- Visual cues for confidence/provenance
- Easy navigation and responsive design

---

## 8. Implementation Steps

1. **Create config files** in `test_config1` and the schema template in `env_templates/`
2. **Add env.testing** with absolute paths
3. **Update env_loader.py** to support absolute path loading and enforce usage in all modules
4. **Refactor extractor.py** for modular source extractors and orchestrator logic
5. **Implement output schema validation**
6. **Scaffold Streamlit UI** for review/edit
7. **Test end-to-end** with dummy and real data

---

## 9. Review & Next Steps

- Review this document and the config/schema samples
- Approve or suggest changes
- Begin implementation step-by-step as outlined above 