# Chemistry Extractor POC: Finalized Design

---

## 1. Folder Structure

```
common_dictionary/
│
├── src/
│   └── domains/
│       └── chemistry/
│           ├── sources/
│           │   ├── base_source.py         # Abstract base for all sources
│           │   ├── pubchem_source.py      # PubChem API extractor
│           │   ├── rdkit_source.py        # RDKit local extractor (to be implemented)
│           │   ├── source_factory.py      # Factory for source instantiation
│           │   └── __init__.py
│           └── extractor.py               # Chemistry-specific orchestrator
│
├── config/
│   └── domains/
│       └── chemistry/
│           └── test_config1/
│               ├── ontology.yaml
│               ├── entity_config.yaml
│               ├── extraction_config.yaml
│               ├── conflict_resolution.yaml
│               ├── relationships_config.yaml
│               ├── validation_config.yaml
│               └── source_mapping.yaml
│
├── env_templates/
│   └── chemistry/
│       └── config_templates/
│           ├── output_schemas.yaml
│           └── ...
│
├── data/
│   └── outputs/
│       └── chemistry/
│           └── test_config1/
│               ├── entities.json
│               └── relationships.json
│
└── documentation/
    └── extractorDesign.md
```

---

## 2. Core Components

### A. Source Extractors (in `sources/`)
- **base_source.py:** Abstract class `BaseSource` with methods:
  - `connect()`
  - `extract_entity(entity_type, query)`
  - `extract_relationship(relationship_type, query)` *(to be added for relationship extraction)*
  - `validate_connection()`
- **pubchem_source.py:** Implements `BaseSource` for PubChem API.
- **rdkit_source.py:** Implements `BaseSource` for local RDKit extraction (to be implemented).
- **source_factory.py:** Factory to instantiate the correct extractor based on config and environment.

### B. Orchestrator (`extractor.py`)
- Loads all configs via `env_loader`.
- For each entity/relationship type in the config:
  - Uses `source_factory` to instantiate the required extractors.
  - Calls `extract_entity` and `extract_relationship` as needed.
  - Merges results, attaches provenance, UUID, and confidence.
  - Handles fallback logic if configured.
- Validates and writes output to the correct directory.

### C. Config Alignment
- All logic, types, and sources are driven by the YAML config files in `test_config1`, which strictly follow the templates in `env_templates/chemistry/config_templates`.

---

## 3. Extraction Models Supported

- **API-Based:** PubChem via `pubchem_source.py`
- **Local Toolkit:** RDKit via `rdkit_source.py`
- **(Optional/Future) ML/LLM-Based:** Can be added as a new extractor class if needed.

---

## 4. Relationship Extraction

- **Config-Driven:** For relationships like `is_subtype_of`, extraction is based on ontology and relationships config.
- **Source-Driven (Future):** If a source provides relationship data, implement `extract_relationship` in the relevant extractor.

---

## 5. Outputs and Output Validation

### A. Output Files
- **entities.json:**
  - List of extracted entities, each with UUID, type, attributes, provenance, confidence, etc.
- **relationships.json:**
  - List of extracted relationships, each with UUID, type, from/to entity references, provenance, confidence, etc.
- **Output Directory:**
  - Set via the environment/config (e.g., `DATA_OUTPUT_DIR` in `env.testing`).

### B. Output Schema Validation
- **Schema:**
  - Defined in `output_schemas.yaml` (from `env_templates/chemistry/config_templates/`).
- **Validation:**
  - Use a library like `jsonschema` to validate each output file (entities and relationships) against the loaded schema.
  - If validation fails, print errors and (optionally) halt the pipeline.
- **Purpose:**
  - Ensures all required fields are present, data types and structure match the template, and outputs are ready for downstream use.

---

## 6. Extensibility

- **Add a new source:** Implement a new class in `sources/`, register it in `source_factory.py`, and update the config.
- **Add a new entity/relationship:** Update the config files; the orchestrator and extractors will pick up the changes automatically if generic.

---

## 7. What's Removed

- **No test_source.py:** All mock/testing logic is removed for production alignment.
- **No example_usage.py:** Not relevant to the orchestrator or extraction pipeline.

---

## 8. Summary Table

| File/Module                | Role/Responsibility                                      |
|----------------------------|---------------------------------------------------------|
| base_source.py             | Abstract base for all extractors                        |
| pubchem_source.py          | PubChem API extractor                                   |
| rdkit_source.py            | RDKit local extractor                                   |
| source_factory.py          | Factory for instantiating extractors by config          |
| extractor.py (chemistry)   | Orchestrates extraction, merging, output, validation    |
| output_schemas.yaml        | Defines required output structure for validation        |
| entities.json, relationships.json | Output files, validated against schema           |
| env_loader.py              | Loads all config paths/environment settings             |
| test_config1/*.yaml        | Drives all logic, types, sources, and mappings          |
| config_templates/*.yaml    | Ensures all configs are aligned and validated           |

---

## 9. Implementation Reference

- This document serves as the implementation reference for the chemistry extractor POC.
- All code and configs should strictly follow this design and the config templates. 