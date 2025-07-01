import os
import logging
from typing import Dict, Any, Optional
from pathlib import Path
from dotenv import load_dotenv

class EnvironmentLoader:
    """Load and manage environment variables for the chemistry domain"""
    
    def __init__(self, domain: str = "chemistry", environment: Optional[str] = None):
        """
        Initialize the environment loader
        
        Args:
            domain: The domain name (e.g., 'chemistry')
            environment: The environment to load (dev, prod, test, staging)
        """
        self.domain = domain
        self.environment = environment or os.getenv('ENVIRONMENT', 'development')
        self.env_vars = {}
        self._load_environment()
        
    def _load_environment(self):
        """Load environment variables from the appropriate .env file"""
        # Get the project root (common_dictionary directory)
        project_root = Path(__file__).parent.parent.parent
        
        # Define the path to env_templates
        env_templates_path = project_root / "env_templates" / self.domain
        
        # Try to load environment-specific file first
        env_file = env_templates_path / f"env.{self.environment}"
        if env_file.exists():
            load_dotenv(env_file, override=True)
            logging.info(f"Loaded environment from {env_file}")
        else:
            # Fall back to template
            template_file = env_templates_path / "env.template"
            if template_file.exists():
                load_dotenv(template_file, override=True)
                logging.info(f"Loaded environment from template {template_file}")
            else:
                logging.warning(f"No environment file found at {env_file} or {template_file}")
        
        self._load_env_vars()
    
    def _load_env_vars(self):
        """Load all relevant environment variables into self.env_vars"""
        # Load all environment variables
        self.env_vars = dict(os.environ)
    
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