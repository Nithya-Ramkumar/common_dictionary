# ===============================================
# CONFIG RECONCILIATION VALIDATION PURPOSES
# ===============================================
# This script performs two critical validations to ensure configuration consistency across domains:
#
# 1. Attribute/Key Presence in entity_config.yaml:
#    - Every attribute listed in 'attributes_to_extract' (in both search-based and key-based mappings)
#      and every key used in key-based mappings must be defined as an attribute for the entity in entity_config.yaml.
#    - Purpose: Ensures that the extraction config is consistent with the entity schema and avoids runtime errors due to undefined attributes.
#
# 2. Key Extractability in Source Mapping:
#    - Every key used in key-based mappings must also be listed as an attribute to extract in the search-based mappings above (or, more generally, must be extractable for the entity).
#    - Purpose: Ensures that the key is actually available for use in key-based lookups, i.e., it is being extracted somewhere in the pipeline and not just referenced.
#
# If either validation fails, the script will raise a clear, actionable error.
# These validations are domain-agnostic and apply to any entity (e.g., Polymer, Compound, etc.).
# ===============================================

import yaml
import os
import sys
import argparse
from config.env_loader import EnvironmentLoader
import logging

logger = logging.getLogger("config_reconciliation")

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
    logger.debug("Reconciling the following YAML config files:")
    logger.debug(f"  Ontology: {ontology_path}")
    logger.debug(f"  Entity config: {entity_config_path}")
    logger.debug(f"  Validation config: {validation_config_path}")
    logger.debug(f"  Source mapping: {source_mapping_path}")
    logger.debug(f"  Conflict resolution: {conflict_resolution_path}")
    logger.debug(f"  Extraction config: {extraction_config_path}")
    logger.debug(f"  Metric units: {metric_units_path}")

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
    logger.debug("\n-------------------------------")
    logger.debug("Ontology Summary")
    logger.debug("-------------------------------")
    logger.debug(f"  ✓ Subdomains found: {', '.join(subdomains) if subdomains else 'None'}")
    logger.debug(f"  ✓ Entity types found: {', '.join(entity_names)}")
    if subclasses:
        logger.debug("  ✓ Subclass relationships:")
        for child, parent in subclasses.items():
            logger.debug(f"    - {child} subclass_of {parent}")
    if subtypes:
        logger.debug("  ✓ Subtypes:")
        for parent, children in subtypes.items():
            logger.debug(f"    - {parent} has subtypes: {', '.join(children)}")
    if tags:
        logger.debug("  ✓ Tags:")
        for entity, taglist in tags.items():
            logger.debug(f"    - {entity}: {', '.join(taglist)}")
    # Print errors or success
    if errors:
        logger.error("\n[Ontology Relationship/Entity Validation Errors]")
        for err in errors:
            logger.error(f"  - {err}")
        raise ValueError("Ontology validation failed. See errors above.")
    else:
        logger.info("[Ontology: all relationships, subclassing, subtypes, and tags are valid]")

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

def extract_attr_names(attr_list):
    """Helper to extract attribute names from a list of strings or dicts."""
    names = set()
    for attr in attr_list:
        if isinstance(attr, dict):
            names.add(attr['name'])
        else:
            names.add(attr)
    return names

def validate_source_mapping_coverage(entity_config, source_mapping):
    """
    Validate that every attribute in entity_config has a corresponding source mapping.
    For the new source_mapping format: check if each attribute appears in any attributes_to_extract
    in search_based_mappings or in any key_based_mappings for that entity.
    """
    logger.debug("\n==============================")
    logger.debug("Source Mapping Coverage Validation")
    logger.debug("==============================")
    
    chem = entity_config.get('chemistry', {})
    entities = {t['name']: t for t in chem.get('types', [])}
    
    missing_mappings = []
    mapped_attrs = []
    
    for entity_name, entity_def in entities.items():
        logger.debug(f"  Entity: {entity_name}")
        attrs = entity_def.get('attributes', [])
        entity_source_mapping = source_mapping.get(entity_name, {})
        # Gather all mapped attributes from search_based_mappings
        search_based = entity_source_mapping.get('search_based_mappings', [])
        search_attrs = set()
        for mapping in search_based:
            search_attrs.update(extract_attr_names(mapping.get('attributes_to_extract', [])))
        # Gather all mapped attributes from key_based_mappings
        key_based = entity_source_mapping.get('key_based_mappings', {})
        key_attrs = set()
        for key, key_map in key_based.items():
            for source in key_map.get('sources', []):
                key_attrs.update(extract_attr_names(source.get('attributes_to_extract', [])))
        all_mapped = search_attrs | key_attrs
        for attr in attrs:
            attr_name = attr['name']
            if attr_name in all_mapped:
                mapped_attrs.append(f"{entity_name}.{attr_name}")
                logger.debug(f"    ✓ {attr_name}: mapped")
            else:
                missing_mappings.append(f"{entity_name}.{attr_name}")
                logger.warning(f"    ✗ {attr_name}: NO SOURCE MAPPING")
    logger.debug("\n==============================")
    logger.debug("Source Mapping Summary")
    logger.debug("==============================")
    logger.debug(f"  ✓ Mapped attributes: {len(mapped_attrs)}")
    logger.debug(f"  ✗ Missing mappings: {len(missing_mappings)}")
    if missing_mappings:
        logger.error("\n[Source Mapping Coverage Errors]")
        logger.error("The following attributes have no source mapping:")
        for missing in missing_mappings:
            logger.error(f"  - {missing}")
        logger.error("\nRecommendations:")
        logger.error("  1. Add source mappings in source_mapping.yaml for these attributes")
        logger.error("  2. Or ensure parent entities have source mappings that can be inherited")
        logger.error("  3. Or remove these attributes if they don't need extraction")
        return False
    else:
        logger.info("\n[Source Mapping Coverage: All attributes have source mappings ✓]")
        return True

def validate_source_mapping_vs_entity_config(entity_config, source_mapping):
    """
    1. Every attribute to extract and every key in source_mapping (both search_based_mappings and key_based_mappings)
       must be defined as an attribute for the entity in entity_config.yaml.
    2. Every key in key_based_mappings must also be listed as an attribute to extract in the search_based_mappings for the same entity.
    """
    errors = []
    chem = entity_config.get('chemistry', {})
    entities = {t['name']: t for t in chem.get('types', [])}

    for entity_name, entity_mapping in source_mapping.items():
        if entity_name not in entities:
            errors.append(f"Entity '{entity_name}' in source_mapping is not defined in entity_config.yaml.")
            continue
        entity_attrs = {attr['name'] for attr in entities[entity_name].get('attributes', [])}
        # 1. Check all attributes to extract in search_based_mappings
        search_based = entity_mapping.get('search_based_mappings', [])
        search_attrs = set()
        for mapping in search_based:
            attrs = mapping.get('attributes_to_extract', [])
            for attr in extract_attr_names(attrs):
                if attr not in entity_attrs:
                    errors.append(f"Attribute '{attr}' in search_based_mappings for entity '{entity_name}' is not defined in entity_config.yaml.")
                search_attrs.add(attr)
        # 2. Check all keys and attributes in key_based_mappings
        key_based = entity_mapping.get('key_based_mappings', {})
        for key, key_mapping in key_based.items():
            # Key must be an attribute in entity_config
            if key not in entity_attrs:
                errors.append(f"Key '{key}' in key_based_mappings for entity '{entity_name}' is not defined as an attribute in entity_config.yaml.")
            # Key must also be extracted in search_based_mappings
            if key not in search_attrs:
                errors.append(f"Key '{key}' in key_based_mappings for entity '{entity_name}' is not listed in attributes_to_extract in any search_based_mappings.")
            # All attributes to extract in key_based sources must be in entity_config
            sources = key_mapping.get('sources', [])
            for source in sources:
                attrs = source.get('attributes_to_extract', [])
                for attr in extract_attr_names(attrs):
                    if attr not in entity_attrs:
                        errors.append(f"Attribute '{attr}' in key_based_mappings for entity '{entity_name}' (key '{key}') is not defined in entity_config.yaml.")
    if errors:
        logger.error("\n[Source Mapping vs Entity Config Validation Errors]")
        for err in errors:
            logger.error(f"  - {err}")
        raise ValueError("Source mapping validation failed. See errors above.")
    else:
        logger.info("[Source mapping: all keys and attributes are valid against entity_config.yaml]")

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

    logger.debug("\n==============================")
    logger.debug("Ontology vs Entity Mapping")
    logger.debug("==============================")
    for ent, meta in entity_map.items():
        logger.debug(f"  ✓ Ontology entity: {ent}")
        if meta['subtypes']:
            logger.debug(f"      - Subclasses: {', '.join(meta['subtypes'])}")
            missing = [s for s in meta['subtypes'] if s not in all_entities]
            if missing:
                logger.warning(f"      ✗ Missing in entity_config: {', '.join(missing)}")
            else:
                logger.debug(f"      - All subclasses present in entity_config: ✓")
        else:
            logger.debug(f"      - Subclasses: None")

    logger.debug("\n==============================")
    logger.debug("Entity Attribute Coverage & Quantitative Field Reconciliation")
    logger.debug("==============================")
    unreconciled_quant = []
    reconciled_quant = []
    for ent, ent_def in all_entities.items():
        logger.debug(f"  Entity: {ent}")
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
            logger.debug(f"    - {name} (type: {typ}) | {'/'.join(status) if status else 'MISSING'}")
    logger.debug("\n==============================")
    logger.debug("Reconciled Quantitative Fields")
    logger.debug("==============================")
    for f in reconciled_quant:
        logger.debug(f"  ✓ {f}")
    logger.debug("\n==============================")
    logger.debug("Unreconciled Quantitative Fields")
    logger.debug("==============================")
    for f in unreconciled_quant:
        logger.warning(f"  ✗ {f}")
    logger.debug("\n==============================")
    logger.debug("Summary")
    logger.debug("==============================")
    logger.debug(f"  ✓ {len(all_entities)} entity types found")
    subclass_count = sum(1 for e in entity_map.values() if e['subclass_of'])
    logger.debug(f"  ✓ {subclass_count} subclasses found and mapped")
    attr_count = sum(len(e.get('attributes', [])) for e in all_entities.values())
    logger.debug(f"  ✓ {attr_count} attributes checked")
    logger.debug(f"  ✓ {len(unreconciled_quant)} unreconciled quantitative fields (see above)")
    
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
    entity_config = load_yaml(entity_config_path)
    source_mapping = load_yaml(source_mapping_path)
    validate_source_mapping_vs_entity_config(entity_config, source_mapping)
    result = exhaustive_config_report()
    return result

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--env', type=str, help='Environment to use (overrides COMMON_DICT_ENV)')
    parser.add_argument('--expected-dir', type=str, help='Expected config directory for verification (overrides EXPECTED_CONFIG_DIR)')
    args, unknown = parser.parse_known_args()
    
    try:
        result = run_config_reconciliation(domain="chemistry", environment=args.env, expected_dir=args.expected_dir)
        
        logger.info("\n==============================")
        logger.info("Final Validation Summary")
        logger.info("==============================")
        logger.info(f"  ✓ Source mapping coverage: {'PASS' if result['source_mapping_coverage_valid'] else 'FAIL'}")
        logger.info(f"  ✓ Unreconciled quantitative fields: {result['unreconciled_quantitative_fields']}")
        
        if result['source_mapping_coverage_valid'] and result['unreconciled_quantitative_fields'] == 0:
            logger.info("\n🎉 All validations passed! Configuration is ready for use.")
            sys.exit(0)
        else:
            logger.warning("\n⚠️  Some validations failed. Please review the issues above.")
            sys.exit(1)
            
    except Exception as e:
        logger.error(f"\n❌ Configuration reconciliation failed: {e}")
        sys.exit(1) 