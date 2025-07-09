# GENERIC CONFIG TEMPLATES FOR DOMAIN EXTRACTION PIPELINES

## Overview (Generic)
These templates are designed to be domain-agnostic and can be used for any structured data extraction pipeline (e.g., chemistry, biology, materials, medicine, etc.).

### How to Use (Generic)
1. **Copy** the templates to your config directory for your domain.
2. **Customize** the fields, entity types, sources, and parameters for your use case.
3. **Align** with your ontology.yaml, which defines the entity type hierarchy, subdomains, tags, and relationships.
4. **Validate** your configuration using the reconciliation script or validation tools.
5. **Run** your extraction pipeline using the customized configs.

### Key Features (Generic)
- **Entity Config:**
  - Supports `subclass_of`, `subtypes`, `tags`, `subdomain`, and `search_term` fields.
  - `search_term` can be used for targeted extraction from external sources (e.g., PubChem, NCBI, etc.).
  - Each attribute can specify type, validation, and provenance requirements.
- **Extraction Config:**
  - Supports any data source: APIs, databases, files, etc.
  - Feature flags and extraction parameters are fully customizable.
  - Per-entity-type extraction limits are supported.
- **Relationships & Conflict Resolution:**
  - Templates support specifying allowed subtypes/tags for relationships and conflict resolution strategies.
- **Runtime Feedback Fields:**
  - Extraction and enrichment outputs will include: `value`, `confidence`, `source`, `version`, and `provenance` for each attribute.

### Example: Using `search_term`
To target a specific category or keyword in an external source, add a `search_term` field to the relevant entity type in your entity config:

```yaml
- name: EntityA
  search_term: "entity_a_keyword"  # Used for targeted extraction from external sources
  ...
```

### Best Practices
- Keep your configs aligned with your ontology and data model.
- Use comments in the templates to guide customization.
- Validate your configs before running extraction.
- Document any domain-specific conventions in your own README.

---
# END GENERIC SECTION

# CHEMISTRY DOMAIN EXAMPLE (Retain below)

## Ontology Alignment
All config templates are designed to align with `ontology.yaml`, which defines the entity type hierarchy (taxonomy), subdomains, tags, and relationships. Templates support `subclass_of`, `subtypes`, and `tags` fields for extensibility and clarity.

## Entity Config Template
- Supports `subclass_of`, `subtypes`, `tags`, and `subdomain` fields.
- Enrichment rules include `input_fields`, `api_url`, `mapping_function`, and `output_fields` for clarity.
- At runtime, each attribute will include feedback fields: `value`, `confidence`, `source`, `version`, and `provenance`.

## Relationships Config Template
- Allows relationships to specify allowed subtypes/tags for source/target entities.
- Aligns with ontology.yaml for consistency.

## Conflict Resolution Template
- Can specify conflict resolution strategies at the subtype level if needed.

## Runtime Feedback Fields
- Extraction and enrichment outputs will include: `value`, `confidence`, `source`, `version`, and `provenance` for each attribute.

## How to Use
- Copy templates to your config directory and customize for your domain.
- Ensure all types, subtypes, and relationships are defined in ontology.yaml.
- Use the reconciliation script to validate consistency.

---
For more details, see the comments in each template file. 