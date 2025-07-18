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

logger = logging.getLogger("extractor")

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

    def flatten_entity_results(self, entity_dict, cid=None, source=None, default_confidence=1.0):
        """
        Converts a dict of attributes (with _provenance) into a list of flat attribute records.
        """
        records = []
        provenance = entity_dict.get("_provenance", {})
        # Use provided cid, or try to extract from entity_dict
        cid_val = cid or entity_dict.get("pubchem_cid") or entity_dict.get("cid")
        for attr, value in entity_dict.items():
            if attr.startswith("_"):
                continue
            records.append({
                "cid": cid_val,
                "attribute": attr,
                "value": value,
                "provenance": provenance,
                "confidence": default_confidence,
                "source": source or provenance.get("source"),
                "timestamp": provenance.get("timestamp", "")
            })
        return records

    def run(self):
        logger.info("Logging is working: Extraction run started.")
        entity_types = self._get_entity_types()
        all_flat_results = []
        errors = []
        for entity_type, entity_schema in entity_types.items():
            mapping = self.source_mapping.get(entity_type, {})
            search_mappings = mapping.get('search_based_mappings', [])
            key_mappings = mapping.get('key_based_mappings', {})
            # --- Search-based extraction ---
            for search_map in search_mappings:
                source_name = search_map['source']
                source_config = next((s for s in self.extraction_config['extraction_sources'] if s['name'] == source_name and s.get('enabled', True)), None)
                if not source_config:
                    logger.warning(f"Source {source_name} not enabled or not found in extraction config.")
                    continue
                source_instance = self.factory.create_source(source_config)
                filters = search_map.get('search_filters', [])
                attributes = search_map.get('attributes_to_extract', [])
                max_results = search_map.get('max_results', 20)
                try:
                    search_results = source_instance.search(
                        entity_type=entity_type,
                        filters=filters,
                        attributes=attributes,
                        max_results=max_results
                    )
                except Exception as e:
                    logger.error(f"Search extraction failed for {source_name}: {e}")
                    errors.append(str(e))
                    search_results = []
                for result in search_results:
                    # Mark missing search-based attrs as 'unavailable'
                    for attr in attributes:
                        if attr not in result:
                            result[attr] = 'unavailable'
                    # --- Key-based extraction ---
                    key_flat_records = []
                    for key, key_map in key_mappings.items():
                        key_value = result.get(key)
                        if not key_value:
                            continue
                        for key_source in key_map.get('sources', []):
                            key_source_name = key_source['name']
                            key_source_config = next((s for s in self.extraction_config['extraction_sources'] if s['name'] == key_source_name and s.get('enabled', True)), None)
                            if not key_source_config:
                                logger.warning(f"Key source {key_source_name} not enabled or not found in extraction config.")
                                continue
                            key_source_instance = self.factory.create_source(key_source_config)
                            key_attrs = key_source.get('attributes_to_extract', [])
                            try:
                                key_result = key_source_instance.extract_by_key(
                                    entity_type=entity_type,
                                    key=key_value,
                                    attributes=key_attrs
                                )
                            except Exception as e:
                                logger.error(f"Key-based extraction failed for {key_source_name}: {e}")
                                errors.append(str(e))
                                key_result = {}
                            # Mark missing key-based attrs as 'unavailable'
                            for attr in key_attrs:
                                attr_name = attr['name'] if isinstance(attr, dict) else attr
                                if attr_name not in key_result:
                                    key_result[attr_name] = 'unavailable'
                            # Flatten key-based results and add to key_flat_records
                            key_flat_records.extend(self.flatten_entity_results(key_result, cid=result.get('pubchem_cid'), source=key_source_name))
                    # Flatten search-based result and add to all_flat_results
                    all_flat_results.extend(self.flatten_entity_results(result, cid=result.get('pubchem_cid'), source=source_name))
                    # Add all key-based flat records
                    all_flat_results.extend(key_flat_records)
        # Output generation (unchanged)
        output_dir = self.env_loader.get('DATA_OUTPUT_DIR', os.getcwd())
        schema_path = self.env_loader.get('OUTPUT_ENTITY_SCHEMA', None)
        if not schema_path or not os.path.exists(schema_path):
            schema_path = os.path.join(os.path.dirname(__file__), '../../../env_templates/chemistry/config_templates/output_entity_schema.template.yaml')
            schema_path = os.path.abspath(schema_path)
        output_gen = OutputGenerator(schema_path)
        summary = {
            'total_results': len(all_flat_results),
            'timestamp': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'errors': errors
        }
        output_gen.generate_outputs(all_flat_results, summary, output_dir, schema_name='default', formats=['json', 'csv', 'tsv', 'html'])
        logger.info("Extraction run complete.")

    def _get_entity_types(self):
        # Assumes chemistry.types structure from entity_config
        chem = self.entity_config.get('chemistry', {})
        return {t['name']: t for t in chem.get('types', [])}

if __name__ == "__main__":
    orchestrator = ExtractorOrchestrator()
    orchestrator.run() 