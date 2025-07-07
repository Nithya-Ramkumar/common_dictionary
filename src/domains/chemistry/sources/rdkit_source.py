from typing import Dict, Any, List
from .base_source import BaseSource
from config.env_loader import EnvironmentLoader

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
except ImportError:
    Chem = None
    Descriptors = None

class RDKitSource(BaseSource):
    """RDKit local cheminformatics extractor"""
    def __init__(self, config, env_loader):
        self.config = config
        self.env_loader = env_loader
        self.connection = config['connection']
        self.schema = config.get('schema', {})
        self.available = Chem is not None

    def connect(self) -> bool:
        """Check if RDKit is available"""
        return self.available

    def extract_entity(self, entity_type, attr_name, endpoint, query):
        # Use schema to determine field logic
        endpoint_info = self.schema.get(endpoint, {})
        fields = endpoint_info.get('fields', [])
        if not self.available:
            return []
        # Find the field config for the requested attribute
        field_config = next((f for f in fields if f['name'].lower() == attr_name.lower()), None)
        if not field_config:
            return []
        # For chemistry, expect a SMILES string in the query
            smiles = query.get('smiles')
        if not smiles:
            return []
                mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return []
        value = None
        if attr_name.lower() == 'canonical_smiles':
            value = Chem.MolToSmiles(mol)
        elif attr_name.lower() == 'molecular_weight':
            value = Descriptors.MolWt(mol)
        # Add more RDKit-driven logic as needed
        if value is not None:
            return [{
                'value': value,
                'provenance': {'source': 'rdkit', 'endpoint': endpoint, 'query': query},
                'confidence': 1.0,
                'source': 'rdkit',
            }]
        return []

    def extract_relationship(self, relationship_type, query):
        # Not implemented for RDKitSource
        return []

    def validate_connection(self) -> bool:
        return self.connect() 