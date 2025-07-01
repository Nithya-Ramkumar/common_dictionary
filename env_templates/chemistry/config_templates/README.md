# Chemistry Config Templates

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