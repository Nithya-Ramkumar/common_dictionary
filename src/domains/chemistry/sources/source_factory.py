from typing import Dict, Any
from ...config.env_loader import EnvironmentLoader
from .base_source import BaseSource
from .test_source import TestSource
from .pubchem_source import PubChemSource

class SourceFactory:
    """Factory class to create appropriate source instances"""
    
    def __init__(self, env_loader: EnvironmentLoader):
        """
        Initialize the source factory
        
        Args:
            env_loader: EnvironmentLoader instance for configuration
        """
        self.env_loader = env_loader
    
    def create_source(self, config: Dict[str, Any]) -> BaseSource:
        """
        Create a source instance based on the config
        
        Args:
            config: Dictionary containing source configuration
                   Must have 'name', 'type', and 'connection' keys
        
        Returns:
            BaseSource: An instance of the appropriate source class
        
        Raises:
            ValueError: If source type is not supported
        """
        source_type = config.get('type')
        name = config.get('name')
        connection = config.get('connection')
        
        if source_type == 'mock':
            return TestSource(name, connection)
        elif source_type == 'pubchem':
            # Check if PubChem is enabled in environment
            if not self.env_loader.get('ENABLE_PUBCHEM', True):
                raise ValueError("PubChem source is disabled in current environment")
            return PubChemSource(name, connection, self.env_loader)
        elif source_type == 'reaxys':
            # Check if Reaxys is enabled in environment
            if not self.env_loader.get('ENABLE_REAXYS', False):
                raise ValueError("Reaxys source is disabled in current environment")
            # TODO: Implement ReaxysSource
            raise ValueError("Reaxys source not yet implemented")
        elif source_type == 'chebi':
            # Check if ChEBI is enabled in environment
            if not self.env_loader.get('ENABLE_CHEBI', False):
                raise ValueError("ChEBI source is disabled in current environment")
            # TODO: Implement ChEBISource
            raise ValueError("ChEBI source not yet implemented")
        
        raise ValueError(f"Unsupported source type: {source_type}")
    
    def create_source_by_name(self, source_name: str) -> BaseSource:
        """
        Create a source instance by name using default configuration
        
        Args:
            source_name: Name of the source to create
        
        Returns:
            BaseSource: An instance of the appropriate source class
        """
        if source_name.lower() == 'pubchem':
            if not self.env_loader.get('ENABLE_PUBCHEM', True):
                raise ValueError("PubChem source is disabled in current environment")
            return PubChemSource(
                name="PubChem",
                connection_string=self.env_loader.get('PUBCHEM_API_BASE_URL'),
                env_loader=self.env_loader
            )
        elif source_name.lower() == 'mock':
            return TestSource("MockSource", "mock://localhost")
        else:
            raise ValueError(f"Unknown source name: {source_name}")
    
    def get_available_sources(self) -> Dict[str, bool]:
        """
        Get a dictionary of available sources and their enabled status
        
        Returns:
            Dict[str, bool]: Dictionary mapping source names to enabled status
        """
        return {
            'pubchem': self.env_loader.get('ENABLE_PUBCHEM', True),
            'reaxys': self.env_loader.get('ENABLE_REAXYS', False),
            'chebi': self.env_loader.get('ENABLE_CHEBI', False),
            'mock': True,  # Mock source is always available
        } 