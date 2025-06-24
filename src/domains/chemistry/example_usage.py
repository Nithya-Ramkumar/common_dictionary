import yaml
import logging
from pathlib import Path
from config.env_loader import EnvironmentLoader
from config.config_loader import ConfigLoader
from sources.source_factory import SourceFactory

def setup_environment(environment: str = "development") -> EnvironmentLoader:
    """Setup environment loader for the specified environment"""
    env_loader = EnvironmentLoader(domain="chemistry", environment=environment)
    logging.info(f"Environment loaded: {environment}")
    return env_loader

def load_sources_with_environment(env_loader: EnvironmentLoader, config_path: str = None):
    """Load and initialize sources using environment configuration"""
    
    # Create source factory with environment loader
    factory = SourceFactory(env_loader)
    
    # If no config path provided, use environment default
    if config_path is None:
        config_path = env_loader.get('EXTRACTION_CONFIG')
    
    # Load extraction config
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    sources = []
    for source_config in config['sources']:
        try:
            source = factory.create_source(source_config)
            if source.connect():
                sources.append(source)
                logging.info(f"Successfully connected to {source.name}")
            else:
                logging.warning(f"Failed to connect to {source_config.get('name', 'unknown')}")
        except Exception as e:
            logging.error(f"Error creating source {source_config.get('name', 'unknown')}: {e}")
    
    return sources

def demonstrate_environment_features(env_loader: EnvironmentLoader):
    """Demonstrate various environment features"""
    
    print("=== Environment Configuration Demo ===")
    print(f"Environment: {env_loader.get('ENVIRONMENT')}")
    print(f"Debug mode: {env_loader.get('DEBUG')}")
    
    # Show available sources
    factory = SourceFactory(env_loader)
    available_sources = factory.get_available_sources()
    print("\nAvailable sources:")
    for source, enabled in available_sources.items():
        status = "✓ Enabled" if enabled else "✗ Disabled"
        print(f"  {source}: {status}")
    
    # Show PubChem settings
    pubchem_settings = env_loader.get_pubchem_settings()
    print(f"\nPubChem settings:")
    for key, value in pubchem_settings.items():
        print(f"  {key}: {value}")
    
    # Show storage paths
    storage_paths = {
        'data_root_dir': env_loader.get('DATA_ROOT_DIR'),
        'chemistry_data_path': env_loader.get('CHEMISTRY_DATA_PATH'),
        'chemistry_raw_data_path': env_loader.get('CHEMISTRY_RAW_DATA_PATH'),
    }
    print(f"\nStorage paths:")
    for key, value in storage_paths.items():
        print(f"  {key}: {value}")
    
    # Show feature flags
    feature_flags = {
        'chemistry_batch_processing': env_loader.get('CHEMISTRY_BATCH_PROCESSING'),
        'chemistry_validation_strict': env_loader.get('CHEMISTRY_VALIDATION_STRICT'),
        'chemistry_synonym_clustering': env_loader.get('CHEMISTRY_SYNONYM_CLUSTERING'),
    }
    print(f"\nFeature flags:")
    for key, value in feature_flags.items():
        status = "✓ Enabled" if value else "✗ Disabled"
        print(f"  {key}: {status}")

def demonstrate_source_usage(env_loader: EnvironmentLoader):
    """Demonstrate source usage with environment configuration"""
    
    print("\n=== Source Usage Demo ===")
    
    # Create source factory
    factory = SourceFactory(env_loader)
    
    # Try to create PubChem source
    try:
        pubchem_source = factory.create_source_by_name('pubchem')
        print(f"✓ Successfully created PubChem source")
        
        # Test connection
        if pubchem_source.connect():
            print("✓ PubChem connection successful")
            
            # Extract some data (example)
            compounds = pubchem_source.extract_entity("compound", {"cid": "2244"})
            print(f"✓ Extracted {len(compounds)} compounds from PubChem")
        else:
            print("✗ PubChem connection failed")
    except Exception as e:
        print(f"✗ Error with PubChem source: {e}")
    
    # Try to create mock source
    try:
        mock_source = factory.create_source_by_name('mock')
        print(f"✓ Successfully created Mock source")
        
        # Extract test data
        compounds = mock_source.extract_entity("TestCompound", {})
        print(f"✓ Extracted {len(compounds)} test compounds from Mock source")
        print(f"  Sample data: {compounds[0] if compounds else 'No data'}")
    except Exception as e:
        print(f"✗ Error with Mock source: {e}")

def demonstrate_config_loader(env_loader: EnvironmentLoader):
    """Demonstrate the enhanced config loader"""
    
    print("\n=== Config Loader Demo ===")
    
    # Create config loader
    config_loader = ConfigLoader(domain="chemistry", environment=env_loader.get('ENVIRONMENT'))
    
    # Load all configs
    try:
        all_configs = config_loader.get_all_configs()
        print("✓ Successfully loaded all configurations")
        
        # Show environment info
        env_info = all_configs.get('environment', {})
        print(f"  Domain: {env_info.get('domain')}")
        print(f"  Environment: {env_info.get('environment')}")
        print(f"  Debug: {env_info.get('debug')}")
        
        # Show feature flags
        feature_flags = env_info.get('feature_flags', {})
        print("  Feature flags:")
        for flag, value in feature_flags.items():
            status = "✓" if value else "✗"
            print(f"    {status} {flag}")
        
        # Show units
        units = env_info.get('units', {})
        print("  Default units:")
        for unit_type, unit in units.items():
            print(f"    {unit_type}: {unit}")
            
    except Exception as e:
        print(f"✗ Error loading configurations: {e}")

# Example usage
if __name__ == "__main__":
    # Setup logging
    logging.basicConfig(level=logging.INFO)
    
    # Demo with development environment
    print("Running demo with development environment...")
    env_loader = setup_environment("development")
    
    # Demonstrate environment features
    demonstrate_environment_features(env_loader)
    
    # Demonstrate source usage
    demonstrate_source_usage(env_loader)
    
    # Demonstrate config loader
    demonstrate_config_loader(env_loader)
    
    print("\n=== Demo Complete ===") 