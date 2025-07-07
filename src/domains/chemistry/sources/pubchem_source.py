from typing import Dict, Any, List
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from config.env_loader import EnvironmentLoader
from .base_source import BaseSource

class PubChemSource(BaseSource):
    """PubChem API data source implementation"""
    
    def __init__(self, config, env_loader):
        self.config = config
        self.env_loader = env_loader
        self.connection = config['connection']
        self.schema = config.get('schema', {})
        self.session = requests.Session()
        self.base_url = self.connection['base_url']
        self.timeout = self.connection.get('timeout', 30)
        self.retries = self.connection.get('retries', 3)
    
    def connect(self) -> bool:
        """Test connection to PubChem API"""
        try:
            response = self.session.get(self.base_url, timeout=self.timeout)
            return response.status_code == 200
        except Exception:
            return False
    
    def extract_entity(self, entity_type, attr_name, endpoint, query):
        """
        Extract entity data from PubChem
        
        Args:
            entity_type: Type of entity to extract (e.g., 'compound', 'substance')
            attr_name: Name of the attribute to extract
            endpoint: Endpoint to use for extraction
            query: Query parameters for the extraction
                  For compounds: {'cid': '2244'} or {'name': 'aspirin'}
        
        Returns:
            List of extracted entities with their properties
        """
        endpoint_info = self.schema.get(endpoint, {})
        id_field = endpoint_info.get('id_field', None)
        fields = endpoint_info.get('fields', [])
        # Find the field config for the requested attribute
        field_config = next((f for f in fields if f['name'].lower() == attr_name.lower()), None)
        if not field_config:
            return []
        # Build the API URL and parameters
        try:
            if endpoint == 'compound':
                # e.g., /compound/name/{name}/JSON
                name = query.get('name')
                if not name:
                    return []
                url = f"{self.base_url}/compound/name/{name}/JSON"
                resp = self.session.get(url, timeout=self.timeout)
                resp.raise_for_status()
                data = resp.json()
                # Try to extract the attribute from the response
                props = data.get('PC_Compounds', [{}])[0]
                # Try to find the field in props
                value = None
                if attr_name.lower() == 'iupacname':
                    # IUPACName is in props['props']
                    for p in props.get('props', []):
                        if p.get('urn', {}).get('label', '').lower() == 'iupac name':
                            value = p.get('value', {}).get('sval')
                            break
                elif attr_name.lower() == 'molecularformula':
                    value = props.get('props', [{}])[0].get('value', {}).get('sval')
                elif attr_name.lower() == 'molecularweight':
                    for p in props.get('props', []):
                        if p.get('urn', {}).get('label', '').lower() == 'molecular weight':
                            value = p.get('value', {}).get('fval')
                            break
                else:
                    # Try to get directly
                    value = props.get(attr_name)
                if value is not None:
                    return [{
                        'value': value,
                        'provenance': {'source': 'pubchem', 'endpoint': endpoint, 'query': query},
                        'confidence': 1.0,
                        'source': 'pubchem',
                    }]
            elif endpoint == 'property':
                # e.g., /compound/cid/{cid}/property/JSON
                cid = query.get('cid')
                if not cid:
                    return []
                url = f"{self.base_url}/compound/cid/{cid}/property/JSON"
                resp = self.session.get(url, timeout=self.timeout)
                resp.raise_for_status()
                data = resp.json()
                props = data.get('PropertyTable', {}).get('Properties', [{}])[0]
                value = props.get(attr_name)
                if value is not None:
                    return [{
                        'value': value,
                        'provenance': {'source': 'pubchem', 'endpoint': endpoint, 'query': query},
                        'confidence': 1.0,
                        'source': 'pubchem',
                    }]
            # Add more endpoints as needed
        except Exception as e:
            # Log error if needed
            return []
            return []
    
    def validate_connection(self) -> bool:
        """Validate that the connection is working"""
        return self.connect() 

    def extract_relationship(self, relationship_type, query):
        # Not implemented for PubChemSource
        return [] 