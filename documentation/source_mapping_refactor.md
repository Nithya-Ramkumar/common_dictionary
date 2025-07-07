# Source Mapping Refactoring: Moving Source References

## Overview
This document describes the refactoring that moved source references from `entity_config.yaml` to `source_mapping.yaml` to improve separation of concerns and maintainability.

## Problem
The original `entity_config.yaml` contained `enrichment_rules` sections that defined source relationships:

```yaml
enrichment_rules:
  - name: pubchem_lookup
    input_fields: [name, formula]
    api_url: "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/..."
    mapping_function: "pubchem_to_compound"
    output_fields: [molecular_mass_value, molecular_mass_unit]
    description: "Look up compound information from PubChem"
    source: "pubchem_api"  # ← This was in the wrong place
    enabled: true
    priority: 1
    parameters: {}
```

## Solution
Moved all source relationships to `source_mapping.yaml` where they belong according to the architecture design:

```yaml
Compound:
  molecular_mass_value:
    - source: reaxys
      endpoint: compound
      priority: 1
    - source: pubchem
      endpoint: property
      priority: 2
  molecular_mass_unit:
    - source: pubchem
      endpoint: property
      priority: 1
    - source: reaxys
      endpoint: compound
      priority: 2
```

## Changes Made

### 1. Updated `entity_config.yaml`
- **Removed**: All `enrichment_rules` sections
- **Kept**: Entity type definitions, attributes, validation rules, and schema information
- **Focus**: What entities and attributes exist, not where to get data from

### 2. Updated `source_mapping.yaml`
- **Added**: Missing entity types (OrganicCompound, Polymer, MetalComplex)
- **Added**: Missing attributes (functional_groups, monomer_units, central_metal)
- **Enhanced**: Proper source mapping with priorities and endpoints for all entity types

### 3. Updated Template Files
- **Updated**: `entity_config.template.yaml` to remove enrichment_rules
- **Updated**: `source_mapping.template.yaml` to include all entity types

## Benefits

### 1. Separation of Concerns
- **Entity Config**: Defines what entities and attributes exist
- **Source Mapping**: Defines where to get data for each attribute

### 2. Single Source of Truth
- All source relationships are now in one place
- Changes to source priorities only require updating `source_mapping.yaml`

### 3. Flexibility
- Same entity can be extracted from multiple sources with different priorities
- Easy to add/remove sources without touching entity definitions

### 4. Maintainability
- Clear separation makes the system easier to understand and modify
- Reduces duplication and potential inconsistencies

## Architecture Alignment

This refactoring aligns with the original design document which states:

> **source_mapping.yaml**: Maps each entity attribute to relevant sources/endpoints, with priority/fallback logic. Decouples mapping/priority logic from extraction source details.

> **entity_config.yaml**: Detailed schema for each entity type: attributes, validation, enrichment, etc. Supports extraction, validation, enrichment, and UI for each entity type.

## Files Modified

1. `common_dictionary/config/domains/chemistry/entity_config.yaml`
2. `common_dictionary/config/domains/chemistry/source_mapping.yaml`
3. `common_dictionary/env_templates/chemistry/config_templates/entity_config.template.yaml`
4. `common_dictionary/env_templates/chemistry/config_templates/source_mapping.template.yaml`

## Validation

The refactored configuration maintains:
- All entity types and attributes from the original config
- All source relationships (moved to appropriate location)
- Proper inheritance structure (subclass_of relationships)
- Template consistency across all files

## Next Steps

1. Update any code that was reading `enrichment_rules` from entity config
2. Ensure the orchestrator uses `source_mapping.yaml` for all source relationships
3. Update documentation to reflect the new structure
4. Test the pipeline with the refactored configuration 