import os
import logging
from typing import Dict, Any, Optional
from pathlib import Path
from dotenv import load_dotenv

# This environment loader picks up the environment setting from the COMMON_DICT_ENV environment variable.
# Set COMMON_DICT_ENV to the desired environment (e.g., 'testing') before running your pipeline.
# The loader will use this value to select the appropriate env.<environment> file from env_templates/<domain>/.
# If COMMON_DICT_ENV is not set, it will default to 'development' and print a warning.
#
# After loading environment variables, the loader will verify that all config file paths (except METRIC_UNITS)
# point to the expected directory (default: 'test_config1'). If not, it will print a warning.
# The expected directory can be set via the EXPECTED_CONFIG_DIR environment variable or the expected_dir argument.
# This ensures that all modules remain agnostic to config file locations and only the loader is responsible for path validation.

class EnvironmentLoader:
    """
    Load and manage environment variables for the chemistry domain.
    - Uses COMMON_DICT_ENV to select the environment (or the environment argument).
    - Verifies that all config file paths (except METRIC_UNITS) point to the expected directory.
    - The expected directory can be set via the EXPECTED_CONFIG_DIR environment variable or the expected_dir argument.
    """
    
    def __init__(self, domain: str = "chemistry", environment: Optional[str] = None, expected_dir: Optional[str] = None):
        """
        Initialize the environment loader
        
        Args:
            domain: The domain name (e.g., 'chemistry')
            environment: The environment to load (overrides COMMON_DICT_ENV if provided)
            expected_dir: The directory name that config paths should contain (default: 'test_config1', can be set via EXPECTED_CONFIG_DIR)
        """
        self.domain = domain
        # Allow expected_dir to be set via environment variable or argument
        self.expected_dir = expected_dir or os.getenv('EXPECTED_CONFIG_DIR', 'test_config1')
        self.environment = (
            environment or os.getenv('COMMON_DICT_ENV') or 'development'
        )
        if not os.getenv('COMMON_DICT_ENV') and not environment:
            print("[env_loader] Warning: COMMON_DICT_ENV not set. Defaulting to 'development'.")
        self.env_vars = {}
        self.debug_flags = {}
        self._load_environment()
        self.verify_config_paths(expected_dir=self.expected_dir)
        
    def _load_environment(self):
        """Load environment variables from the appropriate .env file"""
        # Get the project root (common_dictionary directory)
        project_root = Path(__file__).parent.parent.parent
        
        # Define the path to env_templates
        env_templates_path = project_root / "env_templates" / self.domain
        
        # Try to load environment-specific file first
        env_file = env_templates_path / f"env.{self.environment}"
        print(f"[DEBUG] Environment: {self.environment}")
        print(f"[DEBUG] Attempting to load env file: {env_file}")
        if env_file.exists():
            load_dotenv(env_file, override=True)
            print(f"[DEBUG] Loaded environment from {env_file}")
        else:
            # Fall back to template
            template_file = env_templates_path / "env.template"
            print(f"[DEBUG] Attempting to load template env file: {template_file}")
            if template_file.exists():
                load_dotenv(template_file, override=True)
                print(f"[DEBUG] Loaded environment from template {template_file}")
            else:
                print(f"[DEBUG] No environment file found at {env_file} or {template_file}")
        
        self._load_env_vars()
        # Print out key config paths for verification
        print("[DEBUG] Key config paths loaded:")
        for key in [
            'ONTOLOGY_CONFIG', 'ENTITY_CONFIG', 'EXTRACTION_CONFIG', 'VALIDATION_CONFIG',
            'SOURCE_MAPPING', 'CONFLICT_RESOLUTION', 'RELATIONSHIPS_CONFIG']:
            print(f"  {key}: {os.environ.get(key)}")
    
    def _load_env_vars(self):
        """Load all relevant environment variables into self.env_vars and parse debug flags from environment variables if present"""
        # Load all environment variables
        self.env_vars = dict(os.environ)
        # Load debug flags from environment variables
        self.debug_flags = {
            'pubchem': self.env_vars.get('DEBUG_PUBCHEM', 'false').lower() == 'true',
            'rdkit': self.env_vars.get('DEBUG_RDKIT', 'false').lower() == 'true',
            'reaxys': self.env_vars.get('DEBUG_REAXYS', 'false').lower() == 'true',
        }
    
    def get(self, key: str, default: Any = None) -> Any:
        """Get an environment variable value"""
        return self.env_vars.get(key, default)
    
    def get_all(self) -> Dict[str, Any]:
        """Get all environment variables"""
        return self.env_vars.copy()
    
    def print_config_paths(self):
        """Print all loaded configuration paths"""
        print("All environment variables loaded:")
        for key, value in self.env_vars.items():
            print(f"  {key}: {value}")
    
    def verify_config_paths(self, expected_dir: str = 'test_config1'):
        """
        Verify that all config file paths (except METRIC_UNITS) point to the expected directory.
        Print a warning if any path does not match.
        """
        config_keys = [
            'ONTOLOGY_CONFIG', 'ENTITY_CONFIG', 'EXTRACTION_CONFIG', 'VALIDATION_CONFIG',
            'SOURCE_MAPPING', 'CONFLICT_RESOLUTION', 'RELATIONSHIPS_CONFIG'
        ]
        for key in config_keys:
            path = self.get(key)
            if path and expected_dir not in str(path):
                print(f"[env_loader WARNING] {key} path does not point to {expected_dir}: {path}")
    
    def get_pubchem_settings(self) -> Dict[str, Any]:
        """Get PubChem-specific settings"""
        return {
            'api_base_url': self.get('PUBCHEM_API_BASE_URL'),
            'timeout': self.get('PUBCHEM_TIMEOUT'),
            'max_retries': self.get('PUBCHEM_MAX_RETRIES'),
            'batch_size': self.get('PUBCHEM_BATCH_SIZE'),
        }
    
    def get_database_settings(self) -> Dict[str, Any]:
        """Get database settings"""
        return {
            'host': self.get('DB_HOST'),
            'port': self.get('DB_PORT'),
            'name': self.get('DB_NAME'),
            'user': self.get('DB_USER'),
            'password': self.get('DB_PASSWORD'),
            'pool_size': self.get('DB_POOL_SIZE'),
            'max_overflow': self.get('DB_MAX_OVERFLOW'),
        }
    
    def get_cache_settings(self) -> Dict[str, Any]:
        """Get cache settings"""
        return {
            'type': self.get('CACHE_TYPE'),
            'host': self.get('CACHE_HOST'),
            'port': self.get('CACHE_PORT'),
            'ttl': self.get('CACHE_TTL'),
            'max_size': self.get('CACHE_MAX_SIZE'),
            'enabled': self.get('CACHE_ENABLED'),
        }
    
    def get_logging_settings(self) -> Dict[str, Any]:
        """Get logging settings"""
        return {
            'level': self.get('LOG_LEVEL'),
            'file': self.get('LOG_FILE'),
            'format': self.get('LOG_FORMAT'),
            'rotation': self.get('LOG_ROTATION'),
            'retention': self.get('LOG_RETENTION'),
            'error_log_path': self.get('ERROR_LOG_PATH'),
        }
    
    def is_development(self) -> bool:
        """Check if running in development environment"""
        return self.environment == 'development'
    
    def is_production(self) -> bool:
        """Check if running in production environment"""
        return self.environment == 'production'
    
    def is_testing(self) -> bool:
        """Check if running in testing environment"""
        return self.environment == 'testing'
    
    def is_staging(self) -> bool:
        """Check if running in staging environment"""
        return self.environment == 'staging'
    
    def get_debug_flag(self, source_name: str) -> bool:
        """Get the debug flag for a given source (e.g., 'pubchem') from the YAML config."""
        return bool(self.debug_flags.get(source_name, False)) 