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
            load_dotenv(env_file)
            logging.info(f"Loaded environment from {env_file}")
        else:
            # Fall back to template
            template_file = env_templates_path / "env.template"
            if template_file.exists():
                load_dotenv(template_file)
                logging.info(f"Loaded environment from template {template_file}")
            else:
                logging.warning(f"No environment file found at {env_file} or {template_file}")
        
        # Load all environment variables into self.env_vars
        self._load_env_vars()
    
    def _load_env_vars(self):
        """Load all relevant environment variables into self.env_vars"""
        # Chemistry domain configuration paths
        self.env_vars.update({
            'CONFIG_ROOT': os.getenv('CONFIG_ROOT'),
            'CONFIG_TEST': os.getenv('CONFIG_TEST'),
            'ENTITY_CONFIG': os.getenv('ENTITY_CONFIG'),
            'EXTRACTION_CONFIG': os.getenv('EXTRACTION_CONFIG'),
            'SOURCE_MAPPING': os.getenv('SOURCE_MAPPING'),
            'CONFLICT_RESOLUTION': os.getenv('CONFLICT_RESOLUTION'),
            'VALIDATION_CONFIG': os.getenv('VALIDATION_CONFIG'),
        })
        
        # Chemistry data source API settings
        self.env_vars.update({
            'PUBCHEM_API_BASE_URL': os.getenv('PUBCHEM_API_BASE_URL'),
            'PUBCHEM_TIMEOUT': int(os.getenv('PUBCHEM_TIMEOUT', 30)),
            'PUBCHEM_MAX_RETRIES': int(os.getenv('PUBCHEM_MAX_RETRIES', 3)),
            'PUBCHEM_BATCH_SIZE': int(os.getenv('PUBCHEM_BATCH_SIZE', 50)),
            'REAXYS_API_BASE_URL': os.getenv('REAXYS_API_BASE_URL'),
            'REAXYS_API_KEY': os.getenv('REAXYS_API_KEY'),
            'REAXYS_TIMEOUT': int(os.getenv('REAXYS_TIMEOUT', 30)),
            'REAXYS_MAX_RETRIES': int(os.getenv('REAXYS_MAX_RETRIES', 3)),
            'CHEBI_API_BASE_URL': os.getenv('CHEBI_API_BASE_URL'),
            'CHEBI_TIMEOUT': int(os.getenv('CHEBI_TIMEOUT', 30)),
            'CHEBI_MAX_RETRIES': int(os.getenv('CHEBI_MAX_RETRIES', 3)),
        })
        
        # Rate limiting & retry configuration
        self.env_vars.update({
            'MAX_REQUESTS_PER_MINUTE': int(os.getenv('MAX_REQUESTS_PER_MINUTE', 60)),
            'RATE_LIMIT_WINDOW': int(os.getenv('RATE_LIMIT_WINDOW', 3600)),
            'RETRY_ATTEMPTS': int(os.getenv('RETRY_ATTEMPTS', 3)),
            'RETRY_DELAY': int(os.getenv('RETRY_DELAY', 5)),
            'MAX_CONCURRENT_DOWNLOADS': int(os.getenv('MAX_CONCURRENT_DOWNLOADS', 5)),
        })
        
        # Storage configuration
        self.env_vars.update({
            'DATA_ROOT_DIR': os.getenv('DATA_ROOT_DIR'),
            'RAW_DOCUMENTS_DIR': os.getenv('RAW_DOCUMENTS_DIR'),
            'PROCESSED_DOCUMENTS_DIR': os.getenv('PROCESSED_DOCUMENTS_DIR'),
            'METRICS_DIR': os.getenv('METRICS_DIR'),
            'LOGS_DIR': os.getenv('LOGS_DIR'),
            'CHEMISTRY_DATA_PATH': os.getenv('CHEMISTRY_DATA_PATH'),
            'CHEMISTRY_RAW_DATA_PATH': os.getenv('CHEMISTRY_RAW_DATA_PATH'),
            'CHEMISTRY_PROCESSED_DATA_PATH': os.getenv('CHEMISTRY_PROCESSED_DATA_PATH'),
            'CHEMISTRY_METRICS_PATH': os.getenv('CHEMISTRY_METRICS_PATH'),
        })
        
        # Logging configuration
        self.env_vars.update({
            'LOG_LEVEL': os.getenv('LOG_LEVEL', 'INFO'),
            'LOG_FILE': os.getenv('LOG_FILE'),
            'LOG_FORMAT': os.getenv('LOG_FORMAT'),
            'LOG_ROTATION': os.getenv('LOG_ROTATION'),
            'LOG_RETENTION': os.getenv('LOG_RETENTION'),
            'ERROR_LOG_PATH': os.getenv('ERROR_LOG_PATH'),
        })
        
        # Database configuration
        self.env_vars.update({
            'DB_HOST': os.getenv('DB_HOST', 'localhost'),
            'DB_PORT': int(os.getenv('DB_PORT', 5432)),
            'DB_NAME': os.getenv('DB_NAME'),
            'DB_USER': os.getenv('DB_USER'),
            'DB_PASSWORD': os.getenv('DB_PASSWORD'),
            'DB_POOL_SIZE': int(os.getenv('DB_POOL_SIZE', 5)),
            'DB_MAX_OVERFLOW': int(os.getenv('DB_MAX_OVERFLOW', 10)),
        })
        
        # Cache configuration
        self.env_vars.update({
            'CACHE_TYPE': os.getenv('CACHE_TYPE', 'memory'),
            'CACHE_HOST': os.getenv('CACHE_HOST', 'localhost'),
            'CACHE_PORT': int(os.getenv('CACHE_PORT', 6379)),
            'CACHE_TTL': int(os.getenv('CACHE_TTL', 3600)),
            'CACHE_MAX_SIZE': int(os.getenv('CACHE_MAX_SIZE', 1000)),
            'CACHE_ENABLED': os.getenv('CACHE_ENABLED', 'true').lower() == 'true',
        })
        
        # Security settings
        self.env_vars.update({
            'API_KEY_EXPIRY_DAYS': int(os.getenv('API_KEY_EXPIRY_DAYS', 30)),
            'JWT_SECRET_KEY': os.getenv('JWT_SECRET_KEY'),
            'JWT_ALGORITHM': os.getenv('JWT_ALGORITHM', 'HS256'),
            'TOKEN_EXPIRY_MINUTES': int(os.getenv('TOKEN_EXPIRY_MINUTES', 60)),
            'HASH_ALGORITHM': os.getenv('HASH_ALGORITHM', 'bcrypt'),
        })
        
        # Proxy configuration
        self.env_vars.update({
            'USE_PROXY': os.getenv('USE_PROXY', 'false').lower() == 'true',
            'HTTP_PROXY': os.getenv('HTTP_PROXY'),
            'HTTPS_PROXY': os.getenv('HTTPS_PROXY'),
            'NO_PROXY': os.getenv('NO_PROXY'),
        })
        
        # Error handling
        self.env_vars.update({
            'DETAILED_ERRORS': os.getenv('DETAILED_ERRORS', 'false').lower() == 'true',
            'ALERT_ON_ERROR': os.getenv('ALERT_ON_ERROR', 'true').lower() == 'true',
            'ALERT_EMAIL': os.getenv('ALERT_EMAIL'),
        })
        
        # Feature flags
        self.env_vars.update({
            'ENABLE_RATE_LIMITING': os.getenv('ENABLE_RATE_LIMITING', 'true').lower() == 'true',
            'ENABLE_CACHING': os.getenv('ENABLE_CACHING', 'true').lower() == 'true',
            'ENABLE_MONITORING': os.getenv('ENABLE_MONITORING', 'true').lower() == 'true',
            'ENABLE_ANALYTICS': os.getenv('ENABLE_ANALYTICS', 'true').lower() == 'true',
        })
        
        # Environment settings
        self.env_vars.update({
            'ENVIRONMENT': self.environment,
            'DEBUG': os.getenv('DEBUG', 'false').lower() == 'true',
        })
        
        # Chemistry-specific features
        self.env_vars.update({
            'ENABLE_PUBCHEM': os.getenv('ENABLE_PUBCHEM', 'true').lower() == 'true',
            'ENABLE_REAXYS': os.getenv('ENABLE_REAXYS', 'false').lower() == 'true',
            'ENABLE_CHEBI': os.getenv('ENABLE_CHEBI', 'false').lower() == 'true',
            'CHEMISTRY_BATCH_PROCESSING': os.getenv('CHEMISTRY_BATCH_PROCESSING', 'true').lower() == 'true',
            'CHEMISTRY_VALIDATION_STRICT': os.getenv('CHEMISTRY_VALIDATION_STRICT', 'true').lower() == 'true',
            'CHEMISTRY_SYNONYM_CLUSTERING': os.getenv('CHEMISTRY_SYNONYM_CLUSTERING', 'true').lower() == 'true',
            'DEFAULT_MASS_UNIT': os.getenv('DEFAULT_MASS_UNIT', 'g'),
            'DEFAULT_VOLUME_UNIT': os.getenv('DEFAULT_VOLUME_UNIT', 'mL'),
            'DEFAULT_TEMPERATURE_UNIT': os.getenv('DEFAULT_TEMPERATURE_UNIT', 'K'),
            'DEFAULT_PRESSURE_UNIT': os.getenv('DEFAULT_PRESSURE_UNIT', 'Pa'),
        })
    
    def get(self, key: str, default: Any = None) -> Any:
        """Get an environment variable value"""
        return self.env_vars.get(key, default)
    
    def get_all(self) -> Dict[str, Any]:
        """Get all environment variables"""
        return self.env_vars.copy()
    
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