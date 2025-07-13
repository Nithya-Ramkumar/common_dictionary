from typing import Dict, Any, List
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from config.env_loader import EnvironmentLoader
from .base_source import BaseSource

class PubChemSource(BaseSource):
    """PubChem API data source implementation"""
    
    pubchem_attr_map = {
        'smiles': {
            'pubchem_property': 'SMILES',
            'endpoint': 'property',
            'format': 'json',
            'batchable': True
        },
        'formula': {
            'pubchem_property': 'MolecularFormula',
            'endpoint': 'property',
            'format': 'json',
            'batchable': True
        },
        'iupac_name': {
            'pubchem_property': 'IUPACName',
            'endpoint': 'property',
            'format': 'json',
            'batchable': True
        },
        'molecular_weight': {
            'pubchem_property': 'MolecularWeight',
            'endpoint': 'property',
            'format': 'json',
            'batchable': True
        },
        'synonyms': {
            'pubchem_property': None,
            'endpoint': 'synonyms',
            'format': 'json',
            'batchable': False
        },
        'name': {
            'pubchem_property': 'IUPACName',
            'endpoint': 'property',
            'format': 'json',
            'batchable': True
        },
        'average_molecular_weight': {
            'pubchem_property': 'MolecularWeight',
            'endpoint': 'property',
            'format': 'json',
            'batchable': True
        },
        # Add more as needed ...
    }

    def __init__(self, config, env_loader, debug=False):
        self.config = config
        self.env_loader = env_loader
        self.connection = config['connection']
        self.schema = config.get('schema', {})
        self.session = requests.Session()
        self.base_url = self.connection['base_url']
        self.timeout = self.connection.get('timeout', 30)
        self.retries = self.connection.get('retries', 3)
        self.debug = debug
    
    def connect(self) -> bool:
        """Test connection to PubChem API using a real endpoint."""
        try:
            test_url = f"{self.base_url}/compound/cid/2244/property/MolecularWeight/JSON"
            response = self.session.get(test_url, timeout=self.timeout)
            return response.status_code == 200
        except Exception:
            return False

    def extract_entity(self, entity_type, attr_name, endpoint, query, all_attrs=None, debug_first_cid=None):
        """
        Extract entity data from PubChem using the mapping-driven approach.
        all_attrs: list of all config attribute names to extract (for batch property call)
        debug_first_cid: if True, print the full response and debug info (overrides self.debug if set)
        """
        if all_attrs is None:
            all_attrs = [attr_name]
        cid = query.get('cid')
        if not cid:
            return []
        # Group attributes by endpoint/format/batchable
        batch_groups = {}
        single_groups = []
        for attr in all_attrs:
            config = self.pubchem_attr_map.get(attr.lower())
            if not config:
                continue
            key = (config['endpoint'], config['format'], config['batchable'])
            if config['batchable']:
                batch_groups.setdefault(key, []).append((attr, config['pubchem_property']))
            else:
                single_groups.append((attr, config))
        results = []
        # Batch fetch
        for (endpoint, fmt, _), props in batch_groups.items():
            prop_names = [p for _, p in props]
            prop_list = ','.join(prop_names)
            url = f"{self.base_url}/compound/cid/{cid}/{endpoint}/{prop_list}/{fmt}"
            debug_env = self.env_loader.get('DEBUG_PUBCHEM', False)
            debug = self.debug or debug_env or (debug_first_cid if debug_first_cid is not None else False)
            if debug:
                print(f"[DEBUG] PubChem endpoint: {url}")
                import logging
                logging.info(f"[DEBUG] PubChem endpoint: {url}")
            try:
                resp = self.session.get(url, timeout=self.timeout)
                resp.raise_for_status()
                row = None
                if fmt == 'txt':
                    lines = resp.text.strip().splitlines()
                    if len(lines) < 2:
                        continue
                    headers = lines[0].split('\t')
                    values = lines[1].split('\t')
                    row = dict(zip(headers, values))
                elif fmt == 'csv':
                    import csv
                    import io
                    reader = csv.DictReader(io.StringIO(resp.text))
                    row = next(reader, {})
                elif fmt == 'json':
                    data = resp.json()
                    row = data.get('PropertyTable', {}).get('Properties', [{}])[0]
                else:
                    if debug:
                        print(f"[DEBUG] Unknown format '{fmt}' for PubChem batch property fetch.")
                        logging.info(f"[DEBUG] Unknown format '{fmt}' for PubChem batch property fetch.")
                    continue
                for config_attr, pubchem_prop in props:
                    value = row.get(pubchem_prop) if row else None
                    results.append({
                        'attribute': config_attr,
                        'value': value,
                        'source': 'pubchem',
                        'cid': cid,
                        'provenance': {
                            'source': 'pubchem',
                            'endpoint': endpoint,
                            'query': {'cid': cid}
                        },
                        'confidence': 1.0
                    })
            except Exception as e:
                if debug:
                    print(f"[DEBUG] PubChem batch property fetch failed: {e}")
                    logging.info(f"[DEBUG] PubChem batch property fetch failed: {e}")
                continue
        # Single fetch
        for attr, config in single_groups:
            url = f"{self.base_url}/compound/cid/{cid}/{config['endpoint']}/{config['format']}"
            try:
                resp = self.session.get(url, timeout=self.timeout)
                resp.raise_for_status()
                if config['endpoint'] == 'synonyms' and config['format'] == 'json':
                    data = resp.json()
                    synonyms = data.get('InformationList', {}).get('Information', [{}])[0].get('Synonym', [])
                    results.append({
                        'attribute': attr,
                        'value': synonyms,
                        'provenance': {'source': 'pubchem', 'endpoint': config['endpoint'], 'query': query},
                        'confidence': 1.0,
                        'source': 'pubchem'
                    })
                # Add more endpoint/format handling as needed
            except Exception as e:
                print(f"[ERROR] PubChem extraction failed for CID {cid}, attribute {attr}: {e}")
                continue
        return results
    
    def validate_connection(self) -> bool:
        """Validate that the connection is working"""
        return self.connect() 

    def extract_relationship(self, relationship_type, query):
        # Not implemented for PubChemSource
        return [] 

    def get_cids_by_category(self, term: str, count: int = 20) -> List[str]:
        """Fetch CIDs from PubChem by category search term using ESearch."""
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        params = {
            'db': 'pccompound',
            'term': term,
            'retmax': count,
            'retmode': 'json'
        }
        try:
            resp = self.session.get(url, params=params, timeout=self.timeout)
            resp.raise_for_status()
            data = resp.json()
            return data.get('esearchresult', {}).get('idlist', [])
        except Exception:
            return [] 