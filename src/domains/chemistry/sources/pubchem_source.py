from typing import Dict, Any, List
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from ...config.env_loader import EnvironmentLoader
from .base_source import BaseSource

class PubChemSource(BaseSource):
    """PubChem API data source implementation"""
    
    def __init__(self, name: str, connection_string: str, env_loader: EnvironmentLoader):
        super().__init__(name, connection_string)
        self.env_loader = env_loader
        self.settings = env_loader.get_pubchem_settings()
        self.session = self._setup_session()
        
    def _setup_session(self) -> requests.Session:
        """Setup requests session with retry logic"""
        session = requests.Session()
        retry_strategy = Retry(
            total=self.settings['max_retries'],
            backoff_factor=1,
            status_forcelist=[429, 500, 502, 503, 504]
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        session.mount("https://", adapter)
        session.mount("http://", adapter)
        return session
    
    def connect(self) -> bool:
        """Test connection to PubChem API"""
        try:
            response = self.session.get(
                f"{self.settings['api_base_url']}/ping",
                timeout=self.settings['timeout']
            )
            return response.status_code == 200
        except requests.exceptions.RequestException:
            return False
    
    def extract_entity(self, entity_type: str, query: Dict[str, Any]) -> List[Dict[str, Any]]:
        """
        Extract entity data from PubChem
        
        Args:
            entity_type: Type of entity to extract (e.g., 'compound', 'substance')
            query: Query parameters for the extraction
                  For compounds: {'cid': '2244'} or {'name': 'aspirin'}
        
        Returns:
            List of extracted entities with their properties
        """
        if entity_type.lower() not in ['compound', 'substance']:
            raise ValueError(f"Unsupported entity type for PubChem: {entity_type}")
            
        try:
            # Example: Get compound by CID
            if 'cid' in query:
                response = self.session.get(
                    f"{self.settings['api_base_url']}/compound/cid/{query['cid']}/JSON",
                    timeout=self.settings['timeout']
                )
                if response.status_code == 200:
                    return [self._process_compound_response(response.json())]
            
            # Example: Search compound by name
            elif 'name' in query:
                response = self.session.get(
                    f"{self.settings['api_base_url']}/compound/name/{query['name']}/JSON",
                    timeout=self.settings['timeout']
                )
                if response.status_code == 200:
                    return [self._process_compound_response(response.json())]
            
            return []
            
        except requests.exceptions.RequestException:
            return []
    
    def _process_compound_response(self, response_data: Dict[str, Any]) -> Dict[str, Any]:
        """Process PubChem compound response into standardized format"""
        # Example processing - extend based on your entity schema
        compound = response_data.get('PC_Compounds', [{}])[0]
        return {
            'pubchem_id': str(compound.get('id', {}).get('id', {}).get('cid', '')),
            'molecular_formula': self._get_molecular_formula(compound),
            'molecular_weight': self._get_molecular_weight(compound),
            'iupac_name': self._get_iupac_name(compound),
            'synonyms': self._get_synonyms(compound),
            'source': 'PubChem',
            'source_url': f"https://pubchem.ncbi.nlm.nih.gov/compound/{compound.get('id', {}).get('id', {}).get('cid', '')}"
        }
    
    def _get_molecular_formula(self, compound: Dict[str, Any]) -> str:
        """Extract molecular formula from compound data"""
        props = compound.get('props', [])
        for prop in props:
            if prop.get('urn', {}).get('label') == 'Molecular Formula':
                return prop.get('value', {}).get('sval', '')
        return ''
    
    def _get_molecular_weight(self, compound: Dict[str, Any]) -> float:
        """Extract molecular weight from compound data"""
        props = compound.get('props', [])
        for prop in props:
            if prop.get('urn', {}).get('label') == 'Molecular Weight':
                return float(prop.get('value', {}).get('sval', 0))
        return 0.0
    
    def _get_iupac_name(self, compound: Dict[str, Any]) -> str:
        """Extract IUPAC name from compound data"""
        props = compound.get('props', [])
        for prop in props:
            if prop.get('urn', {}).get('label') == 'IUPAC Name':
                return prop.get('value', {}).get('sval', '')
        return ''
    
    def _get_synonyms(self, compound: Dict[str, Any]) -> List[str]:
        """Extract synonyms from compound data"""
        synonyms = []
        props = compound.get('props', [])
        for prop in props:
            if prop.get('urn', {}).get('label') == 'Synonyms':
                synonyms.extend(prop.get('value', {}).get('sval', []))
        return synonyms
    
    def validate_connection(self) -> bool:
        """Validate that the connection is working"""
        return self.connect() 