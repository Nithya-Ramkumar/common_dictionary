# Extractor Detailed Design

## Overview
This document describes the detailed design of the extractor pipeline for the chemistry domain, focusing on how each configuration file is used, the step-by-step logic, and holistic validation of the system.

---

## 1. Pipeline Overview

The extractor pipeline is modular and config-driven. The main steps are:
1. Load environment and all configs (using the environment loader).
2. Validate config consistency (using config reconciliation).
3. For each entity type to extract:
   - Use source mapping to determine which sources/endpoints to use for each attribute.
   - Use the extractor to call the appropriate source(s) and collect raw data.
   - Use the validator to check data against schema, types, and units.
   - Use conflict resolution to merge/finalize values if multiple sources provide data.
   - Use ontology to ensure entity/relationship structure is valid.
4. Output: Write validated, reconciled entities/relationships to output files.

---

## 2. Config Files and Their Roles

### Ontology Config (`ontology.yaml`)
- Defines entity types, relationships, inheritance, subdomains, and tags.
- Used to determine what entities/relationships to extract and to validate structure.

### Entity Config (`entity_config.yaml`)
- Defines the schema for each entity type: attributes, types, validation rules, enrichment rules.
- Used to know what attributes to extract and how to validate them.

### Extraction Config (`extraction_config.yaml`)
- Lists all available sources, how to connect to them, and (optionally) what data structures/endpoints they provide.
- Used to instantiate source classes and provide connection details.

### Source Mapping Config (`source_mapping.yaml`)
- Maps each entity/attribute to the sources and endpoints to use, with priority/fallback logic.
- Used to determine which source(s) to call for each attribute and in what order.

### Validator Config (`validation_config.yaml`)
- Defines reusable validation rules for attributes.
- Used to validate extracted values.

### Conflict Resolution Config (`conflict_resolution.yaml`)
- Defines how to resolve conflicts when multiple sources provide values for the same attribute.
- Used to merge/finalize values after extraction.

---

## 3. Step-by-Step Extraction Algorithm

1. **Load and Validate All Configs**
   - Use the environment loader to load all config paths.
   - Use config reconciliation to check:
     - All required files are present.
     - All entities/attributes in ontology/entity config are mapped in source mapping and conflict resolution.
     - All sources in source mapping exist in extraction config.

2. **For Each Entity Type (from Ontology/Entity Config)**
   - For each attribute (from entity config):
     - Consult source mapping to get the list of sources/endpoints to use, in order.
     - For each source:
       - Check extraction config for connection details and if enabled.
       - Instantiate the source class using the factory.
       - Call the extraction method (e.g., `extract_entity`) with the right parameters.
       - Collect all values (with provenance, confidence, etc.).

3. **Validate Extracted Data**
   - For each value:
     - Apply validation rules (from entity config and validation config).
     - If invalid, log or discard.

4. **Resolve Conflicts**
   - If multiple values for the same attribute:
     - Apply conflict resolution strategy (from conflict resolution config).
     - If human review is allowed/required, flag for review.

5. **Structure and Output**
   - Ensure output matches ontology/entity config structure.
   - Write output files (entities.json, relationships.json), including provenance, confidence, and review history.

---

## 4. Textual Flowchart

```
Start: Load Environment & Configs
  |
Validate Config Consistency
  |
For Each Entity Type
  |
For Each Attribute
  |
Consult Source Mapping
  |
For Each Source (by priority)
  |
Check Extraction Config & Enabled
  |
Instantiate Source Class (Factory)
  |
Extract Data (raw values, provenance)
  |
Collect All Values
  |
Apply Validation Rules
  |
Apply Conflict Resolution
  |
Flag for Human Review if Needed
  |
Structure Output (entities/relationships)
  |
Write Output Files
  |
End
```

---

## 5. Holistic Config Validation Checklist
- Every entity/attribute in ontology/entity config is mapped in source mapping.
- Every source in source mapping exists and is enabled in extraction config.
- Every attribute has a validation rule in validation config.
- Every attribute has a conflict resolution strategy.
- All schemas and types are consistent across configs.
- No orphaned mappings or missing fields.

---

## 6. Notes
- This design ensures separation of concerns, modularity, traceability, and extensibility.
- Each config file has a single responsibility and is validated both individually and holistically.
- The orchestrator and extractor logic are fully driven by the configs, enabling easy updates and domain extension. 