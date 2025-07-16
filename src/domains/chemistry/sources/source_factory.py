from typing import Dict, Any, Type
from config.env_loader import EnvironmentLoader
# No direct imports of source classes here!

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
        self.env_loader = env_loader
    
    def create_source(self, config: Dict[str, Any]):
        source_type = config.get('type')
        source_class = self._registry.get(source_type)
        if not source_class:
            raise ValueError(f"Unknown source type: {source_type}. Registered types: {list(self._registry.keys())}")
        debug_flag = self.env_loader.get_debug_flag(config.get('name', '').lower())
        return source_class(config, self.env_loader, debug=debug_flag)
    
    def get_available_sources(self) -> Dict[str, bool]:
        return {k: True for k in self._registry.keys()} 