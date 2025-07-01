# Chemistry Common Dictionary Module: Configuration Design

## Overview
This document describes the configuration design for the Chemistry Common Dictionary Module, which provides a modular, extensible, and schema-validated foundation for entity and relationship definitions, extraction, validation, provenance, and human-in-the-loop review.

## Folder Structure
```
common_dictionary/
  config/
    domains/
      chemistry/
        ontology.yaml
        entity_config.yaml
        relationships_config.yaml
        conflict_resolution.yaml
        extraction_config.yaml
        source_mapping.yaml
        validation_config.yaml
  env_templates/
    chemistry/
      config_templates/
        ...
```

## Configuration Files: Roles and Rationale

| File Name                                 | What It Does                                                                 | Why It Exists/Design Rationale                                                                 |
|-------------------------------------------|------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------|
| ontology.yaml                             | Defines the entity type hierarchy (taxonomy), subdomains, tags, and all relationships (with valid source/target types). | Provides a holistic, navigable, and machine-readable view of the domain (ontology + taxonomy). |
| entity_config.yaml                        | Detailed schema for each entity type: attributes, validation, enrichment, etc.; supports subclass_of, subtypes, tags, and subdomain fields. | Supports extraction, validation, enrichment, and UI for each entity type.                      |
| relationships_config.yaml                 | Detailed schema for each relationship type: attributes, validation, etc.; allows relationships to specify allowed subtypes/tags for source/target entities. | Supports extraction, validation, enrichment, and UI for each relationship type.                |
| conflict_resolution.yaml                  | Attribute-level and subtype-level conflict resolution strategies for each entity/attribute.     | Ensures consistent, traceable, and human-in-the-loop conflict handling.                       |
| extraction_config.yaml                    | Lists all extraction sources (APIs, DBs, files), connection details, schema.  | Centralizes extraction source definitions and schemas.                                         |
| source_mapping.yaml                       | Maps each entity attribute to relevant sources/endpoints, with priority/fallback logic. | Decouples mapping/priority logic from extraction source details.                               |
| validation_config.yaml                    | Defines reusable validation rules for attributes, referenced from entity configs. | Ensures consistent, schema-driven validation across the system.                                |

## Configuration Templates
... (rest of the document unchanged) ... 

## Ontology File (`ontology.yaml`)
- **Purpose:** Provides a single source of truth for entity type hierarchy (taxonomy), subdomains, tags, and all relationships.
- **Example:**
```yaml
subdomains:
  - organic
  - inorganic
entity_types:
  - name: Compound
    subtypes: [OrganicCompound, Polymer, MetalComplex]
    tags: [core, chemical]
    subdomain: null
  - name: Polymer
    subclass_of: Compound
    tags: [polymer, macromolecule]
    subdomain: organic
relationships:
  - name: participates_in
    source: Compound
    target: Reaction
    allowed_source_types: [Compound, Polymer, MetalComplex]
    allowed_target_types: [Reaction]
    description: "Compound or subtype participates in a reaction"
```

## New Fields and Features
- **subclass_of, subtypes, tags:** Support for inheritance, extensibility, and flexible grouping.
- **Enrichment rules:** Now include `input_fields`, `api_url`, `mapping_function`, and `output_fields` for clarity and automation.
- **Runtime feedback fields:** Extraction and enrichment outputs include `value`, `confidence`, `source`, `version`, and `provenance` for each attribute.
- **Subtype-aware relationships and conflict resolution:** Relationships and conflict resolution can be defined at the subtype level for fine-grained control.

## Runtime Feedback Loop
- All extracted attributes include feedback fields for traceability and iterative improvement.
- UI and reconciliation scripts are updated to display and validate these fields.

## Example: Entity Config (excerpt)
```yaml
- name: Polymer
  subclass_of: Compound
  tags: [polymer, macromolecule]
  subdomain: organic
  description: "Polymeric chemical compounds"
  attributes:
    - name: monomer_units
      type: list
      required: false
      validation: "list_validation"
      provenance: {}
  enrichment_rules: []
```

## Example: Extraction Output Schema
```json
{
  "entity_type": "Polymer",
  "attributes": {
    "monomer_units": {
      "value": ["styrene", "butadiene"],
      "confidence": 0.95,
      "source": "pubchem",
      "version": "2024-06-10",
      "provenance": {"source": "pubchem", "endpoint": "compound", "extraction_time": "2024-06-10T12:34:56Z"}
    }
  }
}
```

## End-to-End Consistency
- All configs, templates, UI, and scripts are aligned with the ontology and new features.
- The reconciliation script validates ontology alignment, subclassing, enrichment rule structure, and relationship endpoints.
- The Streamlit UI displays and allows editing of all new fields and feedback. 