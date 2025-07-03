from abc import ABC, abstractmethod
from typing import Dict, Any, List
from config.env_loader import EnvironmentLoader

class BaseSource(ABC):
    """Base class for all data sources"""
    
    def __init__(self, name: str, connection_string: str):
        self.name = name
        self.connection_string = connection_string
    
    @abstractmethod
    def connect(self) -> bool:
        """Establish connection to the source"""
        pass
    
    @abstractmethod
    def extract_entity(self, entity_type: str, query: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Extract entity data from the source"""
        pass
    
    @abstractmethod
    def extract_relationship(self, relationship_type: str, query: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Extract relationship data from the source (optional for POC)"""
        pass
    
    @abstractmethod
    def validate_connection(self) -> bool:
        """Validate that the connection is working"""
        pass 