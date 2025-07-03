from pathlib import Path
import os
import logging
from typing import Dict, Any, Optional
import yaml
from dotenv import load_dotenv
from .env_loader import EnvironmentLoader

class ConfigLoader:
    """Loads and manages configuration from environment variables and YAML files"""
    
    def __init__(self, domain: str = "chemistry", environment: Optional[str] = None, env_file: Optional[str] = None):
        """
        Initialize the config loader
        
        Args:
            domain: The domain name (e.g., 'chemistry')
            environment: The environment to load (development, production, testing, staging)
            env_file: Optional path to a specific .env file (for backward compatibility)
        """
        self.domain = domain
        self.environment = environment or os.getenv('ENVIRONMENT', 'development')
        
        # Initialize environment loader
        self.env_loader = EnvironmentLoader(domain=domain, environment=self.environment)
        
        # For backward compatibility, if env_file is provided, load it
        if env_file:
            self._load_env_file(env_file)
        
        self.config_cache: Dict[str, Any] = {}
        self._setup_logging()
    
    def _load_env_file(self, env_file: str) -> None:
        """Load environment variables from a specific .env file (for backward compatibility)"""
        if not os.path.exists(env_file):
            raise FileNotFoundError(f"Environment file not found: {env_file}")
        load_dotenv(env_file)
        logging.info(f"Loaded environment from file: {env_file}")
    
    def _setup_logging(self) -> None:
        """Setup logging based on environment configuration"""
        log_settings = self.env_loader.get_logging_settings()
        log_file = log_settings.get('file')
        log_level = getattr(logging, log_settings.get('level', 'INFO').upper(), logging.INFO)
        log_format = log_settings.get('format', '%(asctime)s - %(name)s - %(levelname)s - %(message)s')

        # Remove all handlers associated with the root logger object
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)

        handlers = []
        if log_file:
            os.makedirs(os.path.dirname(log_file), exist_ok=True)
            handlers.append(logging.FileHandler(log_file))
        handlers.append(logging.StreamHandler())  # Also log to console

        logging.basicConfig(
            level=log_level,
            format=log_format,
            handlers=handlers
        )
    
    def _load_yaml(self, file_path: str) -> Dict[str, Any]:
        """Load and parse a YAML file"""
        if file_path in self.config_cache:
            return self.config_cache[file_path]
            
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Config file not found: {file_path}")
            
        with open(file_path, 'r') as f:
            config = yaml.safe_load(f)
            self.config_cache[file_path] = config
            return config
    
    def get_extraction_config(self) -> Dict[str, Any]:
        """Load the extraction configuration"""
        extraction_config_path = self.env_loader.get('EXTRACTION_CONFIG')
        if not extraction_config_path:
            raise ValueError("EXTRACTION_CONFIG environment variable not set")
        return self._load_yaml(extraction_config_path)
    
    def get_entity_config(self) -> Dict[str, Any]:
        """Load the entity configuration"""
        entity_config_path = self.env_loader.get('ENTITY_CONFIG')
        if not entity_config_path:
            raise ValueError("ENTITY_CONFIG environment variable not set")
        return self._load_yaml(entity_config_path)
    
    def get_source_mapping(self) -> Dict[str, Any]:
        """Load the source mapping configuration"""
        source_mapping_path = self.env_loader.get('SOURCE_MAPPING')
        if not source_mapping_path:
            raise ValueError("SOURCE_MAPPING environment variable not set")
        return self._load_yaml(source_mapping_path)
    
    def get_conflict_resolution(self) -> Dict[str, Any]:
        """Load the conflict resolution configuration"""
        conflict_resolution_path = self.env_loader.get('CONFLICT_RESOLUTION')
        if not conflict_resolution_path:
            raise ValueError("CONFLICT_RESOLUTION environment variable not set")
        return self._load_yaml(conflict_resolution_path)
    
    def get_validation_config(self) -> Dict[str, Any]:
        """Load the validation configuration"""
        validation_config_path = self.env_loader.get('VALIDATION_CONFIG')
        if not validation_config_path:
            raise ValueError("VALIDATION_CONFIG environment variable not set")
        return self._load_yaml(validation_config_path)
    
    def get_pubchem_settings(self) -> Dict[str, Any]:
        """Get PubChem specific settings from environment"""
        return self.env_loader.get_pubchem_settings()
    
    def get_database_settings(self) -> Dict[str, Any]:
        """Get database settings from environment"""
        return self.env_loader.get_database_settings()
    
    def get_cache_settings(self) -> Dict[str, Any]:
        """Get cache settings from environment"""
        return self.env_loader.get_cache_settings()
    
    def get_logging_settings(self) -> Dict[str, Any]:
        """Get logging settings from environment"""
        return self.env_loader.get_logging_settings()
    
    def get_all_configs(self) -> Dict[str, Any]:
        """Load all configuration files"""
        configs = {}
        
        # Load YAML configs
        try:
            configs['entity_config'] = self.get_entity_config()
        except Exception as e:
            logging.warning(f"Failed to load entity_config: {e}")
            configs['entity_config'] = {}
        
        try:
            configs['extraction_config'] = self.get_extraction_config()
        except Exception as e:
            logging.warning(f"Failed to load extraction_config: {e}")
            configs['extraction_config'] = {}
        
        try:
            configs['source_mapping'] = self.get_source_mapping()
        except Exception as e:
            logging.warning(f"Failed to load source_mapping: {e}")
            configs['source_mapping'] = {}
        
        try:
            configs['conflict_resolution'] = self.get_conflict_resolution()
        except Exception as e:
            logging.warning(f"Failed to load conflict_resolution: {e}")
            configs['conflict_resolution'] = {}
        
        try:
            configs['validation_config'] = self.get_validation_config()
        except Exception as e:
            logging.warning(f"Failed to load validation_config: {e}")
            configs['validation_config'] = {}
        
        # Add environment settings
        configs['environment'] = {
            'domain': self.domain,
            'environment': self.environment,
            'debug': self.env_loader.get('DEBUG'),
            'feature_flags': {
                'enable_pubchem': self.env_loader.get('ENABLE_PUBCHEM'),
                'enable_reaxys': self.env_loader.get('ENABLE_REAXYS'),
                'enable_chebi': self.env_loader.get('ENABLE_CHEBI'),
                'chemistry_batch_processing': self.env_loader.get('CHEMISTRY_BATCH_PROCESSING'),
                'chemistry_validation_strict': self.env_loader.get('CHEMISTRY_VALIDATION_STRICT'),
                'chemistry_synonym_clustering': self.env_loader.get('CHEMISTRY_SYNONYM_CLUSTERING'),
            },
            'units': {
                'default_mass_unit': self.env_loader.get('DEFAULT_MASS_UNIT'),
                'default_volume_unit': self.env_loader.get('DEFAULT_VOLUME_UNIT'),
                'default_temperature_unit': self.env_loader.get('DEFAULT_TEMPERATURE_UNIT'),
                'default_pressure_unit': self.env_loader.get('DEFAULT_PRESSURE_UNIT'),
            }
        }
        
        return configs
    
    def get_environment_variable(self, key: str, default: Any = None) -> Any:
        """Get an environment variable value"""
        return self.env_loader.get(key, default)
    
    def get_all_environment_variables(self) -> Dict[str, Any]:
        """Get all environment variables"""
        return self.env_loader.get_all()
    
    def is_development(self) -> bool:
        """Check if running in development environment"""
        return self.env_loader.is_development()
    
    def is_production(self) -> bool:
        """Check if running in production environment"""
        return self.env_loader.is_production()
    
    def is_testing(self) -> bool:
        """Check if running in testing environment"""
        return self.env_loader.is_testing()
    
    def is_staging(self) -> bool:
        """Check if running in staging environment"""
        return self.env_loader.is_staging()
    
    def get_storage_paths(self) -> Dict[str, str]:
        """Get storage paths from environment"""
        return {
            'data_root_dir': self.env_loader.get('DATA_ROOT_DIR'),
            'chemistry_data_path': self.env_loader.get('CHEMISTRY_DATA_PATH'),
            'chemistry_raw_data_path': self.env_loader.get('CHEMISTRY_RAW_DATA_PATH'),
            'chemistry_processed_data_path': self.env_loader.get('CHEMISTRY_PROCESSED_DATA_PATH'),
            'chemistry_metrics_path': self.env_loader.get('CHEMISTRY_METRICS_PATH'),
        }
    
    def get_rate_limiting_settings(self) -> Dict[str, Any]:
        """Get rate limiting settings from environment"""
        return {
            'max_requests_per_minute': self.env_loader.get('MAX_REQUESTS_PER_MINUTE'),
            'rate_limit_window': self.env_loader.get('RATE_LIMIT_WINDOW'),
            'retry_attempts': self.env_loader.get('RETRY_ATTEMPTS'),
            'retry_delay': self.env_loader.get('RETRY_DELAY'),
            'max_concurrent_downloads': self.env_loader.get('MAX_CONCURRENT_DOWNLOADS'),
        } 