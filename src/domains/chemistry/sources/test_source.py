from typing import Dict, Any, List
from .base_source import BaseSource

class TestSource(BaseSource):
    """Mock source for testing purposes"""
    
    def __init__(self, name: str, connection_string: str):
        super().__init__(name, connection_string)
        self.mock_data = {
            "TestCompound": [
                {
                    "test_attribute": "test_value_1",
                    "uuid": "test-compound-001"
                },
                {
                    "test_attribute": "test_value_2",
                    "uuid": "test-compound-002"
                }
            ]
        }
    
    def connect(self) -> bool:
        """Mock connection - always returns True"""
        return True
    
    def extract_entity(self, entity_type: str, query: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Return mock data for the requested entity type"""
        if entity_type in self.mock_data:
            return self.mock_data[entity_type]
        return []
    
    def validate_connection(self) -> bool:
        """Mock validation - always returns True"""
        return True 