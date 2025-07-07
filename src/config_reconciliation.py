# ===============================================
# **CONFIG RECONCILIATION USAGE INSTRUCTIONS**
#
# This script now uses EnvironmentLoader for all config access.
# The loader picks up the environment setting from the COMMON_DICT_ENV environment variable,
# or from the --env command-line argument (which takes precedence).
# The expected config directory for verification can be set via:
#   - the --expected-dir command-line argument (takes precedence)
#   - the EXPECTED_CONFIG_DIR environment variable
#   - or defaults to 'test_config1'
# All config file paths are resolved via env_loader.get(...), which reads from the appropriate env.<environment> file.
#
# Usage examples:
#   export COMMON_DICT_ENV=testing && python3 common_dictionary/src/config_reconciliation.py
#   python3 common_dictionary/src/config_reconciliation.py --env testing
#   python3 common_dictionary/src/config_reconciliation.py --env testing --expected-dir test_config1
#   export EXPECTED_CONFIG_DIR=test_config1 && python3 common_dictionary/src/config_reconciliation.py --env testing
# ===============================================

import yaml
import os
import sys
import argparse
from config.env_loader import EnvironmentLoader

# Option to print all environment variables
PRINT_ENV = os.getenv('PRINT_ENV', '0') == '1' or '--print-env' in sys.argv

# Parse --env and --expected-dir arguments if provided
parser = argparse.ArgumentParser()
parser.add_argument('--env', type=str, help='Environment to use (overrides COMMON_DICT_ENV)')
parser.add_argument('--expected-dir', type=str, help='Expected config directory for verification (overrides EXPECTED_CONFIG_DIR)')
args, unknown = parser.parse_known_args()

# Initialize environment loader (uses --env and --expected-dir if provided, else environment variables or defaults)
env_loader = EnvironmentLoader(domain="chemistry", environment=args.env, expected_dir=args.expected_dir)

# Get config paths from environment only
ontology_path = env_loader.get('ONTOLOGY_CONFIG')
entity_config_path = env_loader.get('ENTITY_CONFIG')
validation_config_path = env_loader.get('VALIDATION_CONFIG')
source_mapping_path = env_loader.get('SOURCE_MAPPING')
conflict_resolution_path = env_loader.get('CONFLICT_RESOLUTION')
extraction_config_path = env_loader.get('EXTRACTION_CONFIG')
metric_units_path = env_loader.get('METRIC_UNITS')

if PRINT_ENV:
    env_loader.print_config_paths()
else:
    print("Reconciling the following YAML config files:")
    print(f"  Ontology: {ontology_path}")
    print(f"  Entity config: {entity_config_path}")
    print(f"  Validation config: {validation_config_path}")
    print(f"  Source mapping: {source_mapping_path}")
    print(f"  Conflict resolution: {conflict_resolution_path}")
    print(f"  Extraction config: {extraction_config_path}")
    print(f"  Metric units: {metric_units_path}")

# Fail fast if any required config path is missing
for var, path in [
    ('ONTOLOGY_CONFIG', ontology_path),
    ('ENTITY_CONFIG', entity_config_path),
    ('VALIDATION_CONFIG', validation_config_path),
    ('SOURCE_MAPPING', source_mapping_path),
    ('CONFLICT_RESOLUTION', conflict_resolution_path),
    ('EXTRACTION_CONFIG', extraction_config_path),
    ('METRIC_UNITS', metric_units_path),
    ('RELATIONSHIPS_CONFIG', env_loader.get('RELATIONSHIPS_CONFIG'))
]:
    if not path:
        raise RuntimeError(f"Missing required config path: {var}")



def validate_ontology_relationships(ontology_path):
    with open(ontology_path, 'r') as f:
        ontology = yaml.safe_load(f)
    entity_names = {e['name'] for e in ontology.get('entity_types', [])}
    subdomains = ontology.get('subdomains', [])
    subclasses = {e['name']: e.get('subclass_of') for e in ontology.get('entity_types', []) if 'subclass_of' in e}
    subtypes = {e['name']: e.get('subtypes', []) for e in ontology.get('entity_types', []) if 'subtypes' in e}
    tags = {e['name']: e.get('tags', []) for e in ontology.get('entity_types', []) if 'tags' in e}
    errors = []
    # Relationship validation
    for rel in ontology.get('relationships', []):
        if rel['source'] not in entity_names:
            errors.append(f"Relationship '{rel['name']}' has invalid source '{rel['source']}' (not an entity type)")
        if rel['target'] not in entity_names:
            errors.append(f"Relationship '{rel['name']}' has invalid target '{rel['target']}' (not an entity type)")
        for src_type in rel.get('allowed_source_types', []):
            if src_type not in entity_names:
                errors.append(f"Relationship '{rel['name']}' allowed_source_type '{src_type}' is not a valid entity type")
        for tgt_type in rel.get('allowed_target_types', []):
            if tgt_type not in entity_names:
                errors.append(f"Relationship '{rel['name']}' allowed_target_type '{tgt_type}' is not a valid entity type")
    # Validate subclass_of, subtypes, tags in entity_types
    for e in ontology.get('entity_types', []):
        if 'subclass_of' in e and e['subclass_of'] not in entity_names:
            errors.append(f"Entity '{e['name']}' has invalid subclass_of '{e['subclass_of']}' (not an entity type)")
        if 'subtypes' in e:
            for subtype in e['subtypes']:
                if subtype not in entity_names:
                    errors.append(f"Entity '{e['name']}' has subtypes entry '{subtype}' that is not a valid entity type")
        if 'tags' in e and not isinstance(e['tags'], list):
            errors.append(f"Entity '{e['name']}' tags field should be a list")
    # Print ontology summary
    print("\n-------------------------------")
    print("Ontology Summary")
    print("-------------------------------")
    print(f"  ✓ Subdomains found: {', '.join(subdomains) if subdomains else 'None'}")
    print(f"  ✓ Entity types found: {', '.join(entity_names)}")
    if subclasses:
        print("  ✓ Subclass relationships:")
        for child, parent in subclasses.items():
            print(f"    - {child} subclass_of {parent}")
    if subtypes:
        print("  ✓ Subtypes:")
        for parent, children in subtypes.items():
            print(f"    - {parent} has subtypes: {', '.join(children)}")
    if tags:
        print("  ✓ Tags:")
        for entity, taglist in tags.items():
            print(f"    - {entity}: {', '.join(taglist)}")
    # Print errors or success
    if errors:
        print("\n[Ontology Relationship/Entity Validation Errors]")
        for err in errors:
            print(f"  - {err}")
        raise ValueError("Ontology validation failed. See errors above.")
    else:
        print("[Ontology: all relationships, subclassing, subtypes, and tags are valid]")

def load_yaml(path):
    with open(path, 'r') as f:
        return yaml.safe_load(f)

def get_entity_and_subclass_map(ontology):
    entity_map = {}
    for e in ontology.get('entity_types', []):
        entity_map[e['name']] = {
            'subclass_of': e.get('subclass_of'),
            'subtypes': e.get('subtypes', []),
            'tags': e.get('tags', []),
            'subdomain': e.get('subdomain'),
        }
    return entity_map

def get_all_entities(entity_config):
    chem = entity_config.get('chemistry', {})
    return {t['name']: t for t in chem.get('types', [])}

def get_validation_rules(validation_config):
    return validation_config.get('entities', {})

def get_source_mapping(source_mapping):
    return source_mapping

def get_metric_units(metric_units):
    # Flatten all allowed units for quick lookup
    allowed_units = set()
    for cat in metric_units.get('base_units', {}).values():
        allowed_units.update(cat)
    for cat in metric_units.get('derived_units', {}).values():
        allowed_units.update(cat)
    return allowed_units

def is_quantitative(attr):
    return attr.get('type') in ('float', 'int') or 'unit' in attr['name']

def validate_source_mapping_coverage(entity_config, source_mapping):
    """
    Validate that every attribute in entity_config has a corresponding source mapping.
    Also checks for inheritance - if a parent entity has source mappings, child entities
    should inherit those mappings for common attributes.
    """
    print("\n==============================")
    print("Source Mapping Coverage Validation")
    print("==============================")
    
    chem = entity_config.get('chemistry', {})
    entities = {t['name']: t for t in chem.get('types', [])}
    
    # Build inheritance map
    inheritance_map = {}
    for entity_name, entity_def in entities.items():
        parent = entity_def.get('subclass_of')
        if parent:
            if parent not in inheritance_map:
                inheritance_map[parent] = []
            inheritance_map[parent].append(entity_name)
    
    missing_mappings = []
    inherited_mappings = []
    direct_mappings = []
    
    for entity_name, entity_def in entities.items():
        print(f"  Entity: {entity_name}")
        attrs = entity_def.get('attributes', [])
        
        for attr in attrs:
            attr_name = attr['name']
            
            # Check if this entity has direct source mapping
            entity_source_mapping = source_mapping.get(entity_name, {})
            has_direct_mapping = attr_name in entity_source_mapping
            
            # Check if parent has source mapping (for inheritance)
            parent = entity_def.get('subclass_of')
            has_inherited_mapping = False
            if parent and parent in source_mapping:
                parent_mapping = source_mapping[parent]
                has_inherited_mapping = attr_name in parent_mapping
            
            # Check if any child entities have this attribute mapped
            children = inheritance_map.get(entity_name, [])
            children_with_mapping = []
            for child in children:
                if child in source_mapping and attr_name in source_mapping[child]:
                    children_with_mapping.append(child)
            
            if has_direct_mapping:
                direct_mappings.append(f"{entity_name}.{attr_name}")
                print(f"    ✓ {attr_name}: direct mapping")
            elif has_inherited_mapping:
                inherited_mappings.append(f"{entity_name}.{attr_name} (inherited from {parent})")
                print(f"    ✓ {attr_name}: inherited from {parent}")
            elif children_with_mapping:
                print(f"    ✓ {attr_name}: mapped in children {', '.join(children_with_mapping)}")
            else:
                missing_mappings.append(f"{entity_name}.{attr_name}")
                print(f"    ✗ {attr_name}: NO SOURCE MAPPING")
    
    print("\n==============================")
    print("Source Mapping Summary")
    print("==============================")
    print(f"  ✓ Direct mappings: {len(direct_mappings)}")
    print(f"  ✓ Inherited mappings: {len(inherited_mappings)}")
    print(f"  ✗ Missing mappings: {len(missing_mappings)}")
    
    if missing_mappings:
        print("\n[Source Mapping Coverage Errors]")
        print("The following attributes have no source mapping (neither direct nor inherited):")
        for missing in missing_mappings:
            print(f"  - {missing}")
        print("\nRecommendations:")
        print("  1. Add source mappings in source_mapping.yaml for these attributes")
        print("  2. Or ensure parent entities have source mappings that can be inherited")
        print("  3. Or remove these attributes if they don't need extraction")
        return False
    else:
        print("\n[Source Mapping Coverage: All attributes have source mappings ✓]")
        return True

def exhaustive_config_report():
    # Use resolved paths
    ontology = load_yaml(ontology_path)
    entity_config = load_yaml(entity_config_path)
    validation_config = load_yaml(validation_config_path)
    source_mapping = load_yaml(source_mapping_path)
    metric_units = load_yaml(metric_units_path)

    entity_map = get_entity_and_subclass_map(ontology)
    all_entities = get_all_entities(entity_config)
    validation_rules = get_validation_rules(validation_config)
    allowed_units = get_metric_units(metric_units)

    print("\n==============================")
    print("Ontology vs Entity Mapping")
    print("==============================")
    for ent, meta in entity_map.items():
        print(f"  ✓ Ontology entity: {ent}")
        if meta['subtypes']:
            print(f"      - Subclasses: {', '.join(meta['subtypes'])}")
            missing = [s for s in meta['subtypes'] if s not in all_entities]
            if missing:
                print(f"      ✗ Missing in entity_config: {', '.join(missing)}")
            else:
                print(f"      - All subclasses present in entity_config: ✓")
        else:
            print(f"      - Subclasses: None")

    print("\n==============================")
    print("Entity Attribute Coverage & Quantitative Field Reconciliation")
    print("==============================")
    unreconciled_quant = []
    reconciled_quant = []
    for ent, ent_def in all_entities.items():
        print(f"  Entity: {ent}")
        attrs = ent_def.get('attributes', [])
        for attr in attrs:
            name = attr['name']
            typ = attr.get('type')
            validation = attr.get('validation')
            is_quant = is_quantitative(attr)
            # Check for validation rule
            has_validation = bool(validation)
            # Check for source mapping
            has_source = name in source_mapping.get(ent, {})
            # Check for unit mapping if quantitative
            has_unit = False
            if is_quant:
                # If it's a unit field, check allowed units
                if 'unit' in name:
                    has_unit = any(u in allowed_units for u in allowed_units)
                else:
                    # For value fields, check if a corresponding unit field exists
                    unit_field = name.replace('_value', '_unit')
                    has_unit = unit_field in [a['name'] for a in attrs]
            status = []
            if has_validation:
                status.append('validation')
            if has_source:
                status.append('source')
            if is_quant:
                if has_unit:
                    status.append('unit')
                    reconciled_quant.append(f"{ent}.{name}")
                else:
                    unreconciled_quant.append(f"{ent}.{name}")
            print(f"    - {name} (type: {typ}) | {'/'.join(status) if status else 'MISSING'}")
    print("\n==============================")
    print("Reconciled Quantitative Fields")
    print("==============================")
    for f in reconciled_quant:
        print(f"  ✓ {f}")
    print("\n==============================")
    print("Unreconciled Quantitative Fields")
    print("==============================")
    for f in unreconciled_quant:
        print(f"  ✗ {f}")
    print("\n==============================")
    print("Summary")
    print("==============================")
    print(f"  ✓ {len(all_entities)} entity types found")
    subclass_count = sum(1 for e in entity_map.values() if e['subclass_of'])
    print(f"  ✓ {subclass_count} subclasses found and mapped")
    attr_count = sum(len(e.get('attributes', [])) for e in all_entities.values())
    print(f"  ✓ {attr_count} attributes checked")
    print(f"  ✓ {len(unreconciled_quant)} unreconciled quantitative fields (see above)")
    
    # Validate source mapping coverage
    source_mapping_coverage_valid = validate_source_mapping_coverage(entity_config, source_mapping)
    
    return {
        'source_mapping_coverage_valid': source_mapping_coverage_valid,
        'unreconciled_quantitative_fields': len(unreconciled_quant)
    }

def run_config_reconciliation(domain, environment=None, expected_dir=None):
    env_loader = EnvironmentLoader(domain=domain, environment=environment, expected_dir=expected_dir)
    ontology_path = env_loader.get('ONTOLOGY_CONFIG')
    entity_config_path = env_loader.get('ENTITY_CONFIG')
    validation_config_path = env_loader.get('VALIDATION_CONFIG')
    source_mapping_path = env_loader.get('SOURCE_MAPPING')
    conflict_resolution_path = env_loader.get('CONFLICT_RESOLUTION')
    extraction_config_path = env_loader.get('EXTRACTION_CONFIG')
    metric_units_path = env_loader.get('METRIC_UNITS')
    relationships_config_path = env_loader.get('RELATIONSHIPS_CONFIG')
    # Fail fast if any required config path is missing
    for var, path in [
        ('ONTOLOGY_CONFIG', ontology_path),
        ('ENTITY_CONFIG', entity_config_path),
        ('VALIDATION_CONFIG', validation_config_path),
        ('SOURCE_MAPPING', source_mapping_path),
        ('CONFLICT_RESOLUTION', conflict_resolution_path),
        ('EXTRACTION_CONFIG', extraction_config_path),
        ('METRIC_UNITS', metric_units_path),
        ('RELATIONSHIPS_CONFIG', relationships_config_path)
    ]:
        if not path:
            raise RuntimeError(f"Missing required config path: {var}")
    validate_ontology_relationships(ontology_path)
    result = exhaustive_config_report()
    return result

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--env', type=str, help='Environment to use (overrides COMMON_DICT_ENV)')
    parser.add_argument('--expected-dir', type=str, help='Expected config directory for verification (overrides EXPECTED_CONFIG_DIR)')
    args, unknown = parser.parse_known_args()
    
    try:
        result = run_config_reconciliation(domain="chemistry", environment=args.env, expected_dir=args.expected_dir)
        
        print("\n==============================")
        print("Final Validation Summary")
        print("==============================")
        print(f"  ✓ Source mapping coverage: {'PASS' if result['source_mapping_coverage_valid'] else 'FAIL'}")
        print(f"  ✓ Unreconciled quantitative fields: {result['unreconciled_quantitative_fields']}")
        
        if result['source_mapping_coverage_valid'] and result['unreconciled_quantitative_fields'] == 0:
            print("\n🎉 All validations passed! Configuration is ready for use.")
            sys.exit(0)
        else:
            print("\n⚠️  Some validations failed. Please review the issues above.")
            sys.exit(1)
            
    except Exception as e:
        print(f"\n❌ Configuration reconciliation failed: {e}")
        sys.exit(1) 