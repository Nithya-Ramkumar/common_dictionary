from abc import ABC, abstractmethod
from typing import Dict, Any, List, Optional
from config.env_loader import EnvironmentLoader

class BaseSource(ABC):
    """
    Base class for all data sources.
    Defines the interface for search-based and key-based extraction.
    Subclasses must implement the search and extract_by_key methods.
    """
    # Mapping from entity attribute names to source-specific field names
    attribute_map: Dict[str, str] = {}

    def __init__(self, config: Dict[str, Any], env_loader: Optional[EnvironmentLoader] = None, debug: bool = False):
        self.config = config
        self.env_loader = env_loader
        self.debug = debug
        self.name = config.get('name', '')
        self.connection_string = config.get('connection', {}).get('base_url', '')

    @abstractmethod
    def search(self, entity_type: str, filters: List[Dict[str, Any]], attributes: List[str], max_results: int) -> List[Dict[str, Any]]:
        """
        Perform search-based extraction.
        Args:
            entity_type: Name of the entity type (e.g., 'Polymer')
            filters: List of filter dicts (from source mapping)
            attributes: List of attributes to extract
            max_results: Maximum number of results to return
        Returns:
            List of dicts, each with all requested attributes (missing as 'unavailable')
        """
        pass

    @abstractmethod
    def extract_by_key(self, entity_type: str, key: Any, attributes: List[str]) -> Dict[str, Any]:
        """
        Perform key-based extraction for a single key (or a batch if supported).
        Args:
            entity_type: Name of the entity type
            key: The key value (e.g., SMILES string)
            attributes: List of attributes to extract
        Returns:
            Dict with all requested attributes (missing as 'unavailable')
        """
        pass 