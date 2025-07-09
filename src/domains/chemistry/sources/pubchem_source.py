from typing import Dict, Any, List
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from config.env_loader import EnvironmentLoader
from .base_source import BaseSource

class PubChemSource(BaseSource):
    """PubChem API data source implementation"""
    
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
            # Use a known CID and property for a lightweight check
            test_url = f"{self.base_url}/compound/cid/2244/property/MolecularWeight/JSON"
            response = self.session.get(test_url, timeout=self.timeout)
            return response.status_code == 200
        except Exception:
            return False
    
    def extract_entity(self, entity_type, attr_name, endpoint, query, all_attrs=None, debug_first_cid=None):
        """
        Extract entity data from PubChem using the /property/ endpoint for all mapped attributes.
        all_attrs: list of all config attribute names to extract (for batch property call)
        debug_first_cid: if True, print the full response and debug info (overrides self.debug if set)
        """
        # Attribute name mapping for PubChem
        pubchem_attr_map = {
            'average_molecular_weight': 'MolecularWeight',
            'molecular_weight': 'MolecularWeight',
            'iupac_name': 'IUPACName',
            'iupacname': 'IUPACName',
            'formula': 'MolecularFormula',
            'smiles': 'CanonicalSMILES',
            'name': 'IUPACName',
            # Add more mappings as needed
        }
        # If all_attrs is not provided, just extract the single attr_name
        if all_attrs is None:
            all_attrs = [attr_name]
        pubchem_props = [pubchem_attr_map.get(a.lower(), a) for a in all_attrs]
        cid = query.get('cid')
        if not cid:
            return []
        prop_list = ','.join(pubchem_props)
        url = f"{self.base_url}/compound/cid/{cid}/property/{prop_list}/JSON"
        try:
            resp = self.session.get(url, timeout=self.timeout)
            resp.raise_for_status()
            data = resp.json()
            debug = self.debug if debug_first_cid is None else debug_first_cid
            if debug:
                print(f"[DEBUG] PubChem response for CID {cid}:\n{data}")
            props = data.get('PropertyTable', {}).get('Properties', [{}])[0]
            results = []
            for config_attr in all_attrs:
                pubchem_attr = pubchem_attr_map.get(config_attr.lower(), config_attr)
                value = props.get(pubchem_attr)
                if value is not None:
                    results.append({
                        'attribute': config_attr,
                        'value': value,
                        'provenance': {'source': 'pubchem', 'endpoint': 'property', 'query': query},
                        'confidence': 1.0,
                        'source': 'pubchem',
                    })
                else:
                    if debug:
                        print(f"[DEBUG] Attribute '{config_attr}' (PubChem: '{pubchem_attr}') missing for CID {cid}")
            return results
        except Exception as e:
            print(f"[ERROR] PubChem extraction failed for CID {cid}, properties {prop_list}: {e}")
            return []
    
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