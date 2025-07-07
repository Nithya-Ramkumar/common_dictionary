from typing import Dict, Any, Type
from config.env_loader import EnvironmentLoader
from domains.chemistry.sources.pubchem_source import PubChemSource
from domains.chemistry.sources.rdkit_source import RDKitSource
from domains.chemistry.sources.base_source import BaseSource

class SourceFactory:
    """
    Generic, extensible factory for source extractors.
    Uses a registry/plugin pattern: each source class registers itself with a unique type string.
    The factory instantiates the correct class based on the 'type' field in the config.
    """
    _registry: Dict[str, Type] = {}

    @classmethod
    def register_source_type(cls, source_type: str, source_class: Type):
        cls._registry[source_type] = source_class
    
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
        source_class = self._registry.get(source_type)
        if not source_class:
            raise ValueError(f"Unknown source type: {source_type}. Registered types: {list(self._registry.keys())}")
        return source_class(config, self.env_loader)
    
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
        elif source_name.lower() == 'rdkit':
            return RDKitSource("RDKit", "local")
        else:
            raise ValueError(f"Unknown source name: {source_name}")
    
    def get_available_sources(self) -> Dict[str, bool]:
        """
        Get a dictionary of available sources and their enabled status
        
        Returns:
            Dict[str, bool]: Dictionary mapping source names to enabled status
        """
        return {k: True for k in self._registry.keys()}

# Example: Registering source types (to be done in each source module)
# from .pubchem_source import PubChemSource
# SourceFactory.register_source_type("pubchem", PubChemSource)
# from .rdkit_source import RDKitSource
# SourceFactory.register_source_type("rdkit", RDKitSource)
# ...register more as needed... 