import os
import json
import uuid
from typing import List, Dict, Any
import yaml
import jsonschema
import argparse
from domains.chemistry.sources.source_factory import SourceFactory
from config.env_loader import EnvironmentLoader
from config.config_loader import ConfigLoader
import logging
# Ensure all source types are registered before use
import domains.chemistry.sources
from domains.chemistry.output_generator import OutputGenerator
import datetime

# ===============================================
# **CHEMISTRY EXTRACTOR USAGE INSTRUCTIONS**
#
# This script uses EnvironmentLoader for all config access.
# The loader picks up the environment setting from the COMMON_DICT_ENV environment variable,
# or from the --env command-line argument (which takes precedence).
#
# Usage examples:
#   export COMMON_DICT_ENV=testing && python3 common_dictionary/src/domains/chemistry/extractor.py
#   python3 common_dictionary/src/domains/chemistry/extractor.py --env testing
#
# Note: The expected_dir/EXPECTED_CONFIG_DIR logic is only used in config_reconciliation.py for reconciliation diagnostics.
#       It is NOT used in the extractor.
# ===============================================

class ExtractorOrchestrator:
    def __init__(self, domain="chemistry", environment=None):
        self.config_loader = ConfigLoader(domain=domain, environment=environment)
        self.env_loader = self.config_loader.env_loader
        self.extraction_config = self.config_loader.get_extraction_config()
        self.entity_config = self.config_loader.get_entity_config()
        self.source_mapping = self.config_loader.get_source_mapping()
        self.validation_config = self.config_loader.get_validation_config()
        self.conflict_resolution = self.config_loader.get_conflict_resolution()
        self.ontology = self.config_loader.get_environment_variable('ONTOLOGY_CONFIG')
        self.factory = SourceFactory(self.env_loader)
        self.entities = []
        self.relationships = []

    def run(self):
        logging.info("Logging is working: Extraction run started.")
        # Check PubChem API reachability
        pubchem_config = next((s for s in self.extraction_config['extraction_sources'] if s['name'] == 'pubchem' and s.get('enabled', True)), None)
        pubchem_reachable = False
        if pubchem_config:
            from domains.chemistry.sources.pubchem_source import PubChemSource
            pubchem_source = PubChemSource(pubchem_config, self.env_loader)
            pubchem_reachable = pubchem_source.connect()
            if pubchem_reachable:
                logging.info("PubChem API is reachable.")
            else:
                logging.warning("PubChem API is NOT reachable!")
        else:
            logging.warning("PubChem source config not found or not enabled.")

        # Utility functions for fetching CIDs and names
        def fetch_random_cids(count=20, type_='compound'):
            import requests
            try:
                url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/JSON?list_return=listkey&listkey_count={count}"
                resp = requests.get(url, timeout=10)
                resp.raise_for_status()
                data = resp.json()
                cids = data.get('IdentifierList', {}).get('CID', [])
                return cids[:count]
            except Exception as e:
                logging.warning(f"Failed to fetch random CIDs from PubChem: {e}")
                return []

        def fetch_names_for_cids(cids):
            import requests
            names = []
            for cid in cids:
                try:
                    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IUPACName/JSON"
                    resp = requests.get(url, timeout=10)
                    resp.raise_for_status()
                    data = resp.json()
                    props = data.get('PropertyTable', {}).get('Properties', [{}])[0]
                    name = props.get('IUPACName')
                    if name:
                        names.append(name)
                except Exception:
                    continue
            return names

        # Get number of entities per type from config (default 20)
        extraction_params = self.extraction_config.get('extraction_parameters', {})
        max_entities_param = extraction_params.get('max_entities_per_type', 20)
        # max_entities_param can be an int (global) or a dict (per-entity-type)
        def get_limit_for_entity(entity_type):
            if isinstance(max_entities_param, dict):
                return max_entities_param.get(entity_type, max_entities_param.get('default', 20))
            return max_entities_param
        # Get entity types
        entity_types = self._get_entity_types()
        # Track counts by entity and by source
        entity_type_counts = {}
        source_counts = {}
        field_errors = {}
        exceptions = []
        for entity_type, entity_schema in entity_types.items():
            limit = get_limit_for_entity(entity_type)
            logging.info(f"Processing entity type: {entity_type} (limit: {limit})")
            attributes = entity_schema.get('attributes', [])
            search_term = entity_schema.get('search_term', entity_type)
            if entity_type == 'Polymer':
                if pubchem_config:
                    from domains.chemistry.sources.pubchem_source import PubChemSource
                    pubchem_source = PubChemSource(pubchem_config, self.env_loader)
                    cids = pubchem_source.get_cids_by_category(search_term, limit)
                    queries = [ {'cid': cid} for cid in cids ]
                else:
                    queries = []
            elif entity_type == 'Compound':
                cids = fetch_random_cids(limit, type_=entity_type.lower())
                names = fetch_names_for_cids(cids)
                queries = []
                for name, cid in zip(names, cids):
                    queries.append({'name': name, 'cid': cid})
                if not queries:
                    queries = [{'cid': cid} for cid in cids]
                queries = queries[:limit]
            else:
                queries = [{}]
            for query in queries:
                for attr in attributes:
                    attr_name = attr['name']
                    logging.info(f"  Attribute: {attr_name}")
                    mapping = self.source_mapping.get(entity_type, {}).get(attr_name, [])
                    values = []
                    error_for_field = None
                    for source_map in sorted(mapping, key=lambda x: x.get('priority', 1)):
                        if not source_map.get('enabled', True):
                            continue
                        source_name = source_map['source']
                        endpoint = source_map['endpoint']
                        source_config = next((s for s in self.extraction_config['extraction_sources'] if s['name'] == source_name and s.get('enabled', True)), None)
                        if not source_config:
                            logging.warning(f"    Source {source_name} not enabled or not found in extraction config.")
                            continue
                        source_instance = self.factory.create_source(source_config)
                        logging.info(f"    Querying source: {source_name}, endpoint: {endpoint}, query: {query}")
                        try:
                            extracted = source_instance.extract_entity(entity_type, attr_name, endpoint, query)
                            logging.info(f"    Extracted: {extracted}")
                        except Exception as e:
                            error_for_field = str(e)
                            exceptions.append(f"{entity_type}.{attr_name}: {e}")
                            extracted = []
                        if extracted:
                            values.extend(extracted)
                            entity_type_counts.setdefault(entity_type, 0)
                            entity_type_counts[entity_type] += 1
                            source_counts.setdefault(source_name, 0)
                            source_counts[source_name] += 1
                            if not source_map.get('fallback', False):
                                break
                    validated = values
                    final_value = validated[0] if validated else None
                    if final_value:
                        entity_record = {
                            'entity_type': entity_type,
                            'attribute': attr_name,
                            'value': final_value.get('value'),
                            'provenance': final_value.get('provenance'),
                            'confidence': final_value.get('confidence'),
                            'source': final_value.get('source'),
                            'timestamp': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                        }
                        prov = final_value.get('provenance', {})
                        query_cid = prov.get('query', {}).get('cid')
                        if query_cid:
                            entity_record['cid'] = query_cid
                        self.entities.append(entity_record)
                    else:
                        # Track field errors if extraction failed
                        prov = query.get('cid', None)
                        field_errors.setdefault(attr_name, []).append(prov)
                        if error_for_field:
                            field_errors[attr_name][-1] = error_for_field
        # Save summary counts and errors for output
        unique_cids = set(ent['cid'] for ent in self.entities if 'cid' in ent)
        # Get requested count from extraction config if available
        extraction_params = self.extraction_config.get('extraction_parameters', {})
        max_entities_param = extraction_params.get('max_entities_per_type', 20)
        if isinstance(max_entities_param, dict):
            polymers_requested = max_entities_param.get('Polymer', max_entities_param.get('default', 20))
        else:
            polymers_requested = max_entities_param
        summary = {
            'polymers_requested': polymers_requested,
            'polymers_extracted': len(unique_cids),
            'timestamp': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'total_entities': sum(entity_type_counts.values()),
            'entity_type_counts': entity_type_counts,
            'source_counts': source_counts,
            'field_errors': field_errors,
            'exceptions': exceptions
        }
        # Determine schema path (prefer config dir, fallback to template)
        schema_path = os.path.join(self.env_loader.get('CONFIG_DIR', ''), 'output_entity_schema.yaml')
        if not os.path.exists(schema_path):
            # Fallback to template
            schema_path = os.path.join(os.path.dirname(__file__), '../../../env_templates/chemistry/config_templates/output_entity_schema.template.yaml')
            schema_path = os.path.abspath(schema_path)
        output_dir = self.env_loader.get('DATA_OUTPUT_DIR', os.getcwd())
        output_gen = OutputGenerator(schema_path)
        output_gen.generate_outputs(self.entities, summary, output_dir, schema_name='default', formats=['json', 'csv', 'tsv', 'html'])
        logging.info("Extraction run complete.")

    def _get_entity_types(self):
        # Assumes chemistry.types structure from entity_config
        chem = self.entity_config.get('chemistry', {})
        return {t['name']: t for t in chem.get('types', [])}

if __name__ == "__main__":
    orchestrator = ExtractorOrchestrator()
    orchestrator.run() 