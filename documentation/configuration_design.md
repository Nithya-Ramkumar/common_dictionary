# Chemistry Common Dictionary Module: Configuration Design

## Overview
This document describes the configuration design for the Chemistry Common Dictionary Module, which provides a modular, extensible, and schema-validated foundation for entity and relationship definitions, extraction, validation, provenance, and human-in-the-loop review.

## Folder Structure
```
common_dictionary/
  env_templates/
    chemistry/
      config_templates/
        entity_config.template.yaml
        extraction_config.template.yaml
        source_mapping.template.yaml
        validation_config.template.yaml
        conflict_resolution.template.yaml
        relationships_config.template.yaml
```

## Configuration Templates

### 1. `entity_config.template.yaml`
- **Purpose:** Defines all entity types, their attributes, validation rules, provenance tracking, and enrichment rules.
- **Key Features:**
  - Domain and subdomain tagging
  - Attribute-level schema and validation
  - Provenance fields for each attribute
  - Enrichment rules for external data sources
  - Extensible for subdomains and overrides

### 2. `extraction_config.template.yaml`
- **Purpose:** Lists all extraction sources (APIs, databases, files), their connection details, and response schemas.
- **Key Features:**
  - Source type (API, database, file)
  - Connection/authentication details
  - Endpoint/table/file schema
  - No mapping or priority logic (see source_mapping)

### 3. `source_mapping.template.yaml`
- **Purpose:** Maps each entity attribute to all relevant sources and endpoints, with priority and fallback logic.
- **Key Features:**
  - Attribute-level source mapping
  - Priority and fallback configuration
  - Endpoint and parameter specification
  - No connection details (see extraction_config)

### 4. `validation_config.template.yaml`
- **Purpose:** Defines reusable validation rules for attributes, referenced from entity configs.
- **Key Features:**
  - Patterns, types, required flags, ranges, enums
  - Unit validation (linked to `metric_units.yaml`)
  - Custom logic for domain-specific validation

### 5. `conflict_resolution.template.yaml`
- **Purpose:** Specifies conflict resolution strategies for each entity attribute when multiple sources disagree.
- **Key Features:**
  - Attribute-level strategy (e.g., source priority, average, union)
  - Human-in-the-loop escalation options
  - Notes for custom logic

### 6. `relationships_config.template.yaml`
- **Purpose:** Defines relationship types between entities, including directionality, allowed types, and validation.
- **Key Features:**
  - Relationship type and description
  - Source/target entity types
  - Directionality and cardinality
  - Validation rules

## Usage Workflow
1. **Copy templates** to your domain directory and customize as needed.
2. **Define entities and attributes** in `entity_config.yaml`.
3. **List extraction sources** in `extraction_config.yaml`.
4. **Map sources to attributes** in `source_mapping.yaml`.
5. **Specify validation rules** in `validation_config.yaml`.
6. **Configure conflict resolution** in `conflict_resolution.yaml`.
7. **Define relationships** in `relationships_config.yaml`.
8. **Test all configs** with `config_reconciliation.py` before deployment.
9. **Iterate** with human-in-the-loop review and export validated outputs.

## Provenance and Human-in-the-Loop
- All sources for each attribute are tracked and displayed in the UI.
- Human reviewers can resolve conflicts, validate data, and export results.
- Provenance is maintained for traceability and iterative improvement.

## Extensibility
- Add new domains by copying and customizing templates.
- Add subdomains or override configs as needed.
- Extend validation and enrichment rules for new entity types.

---
**End of configuration design documentation.** 