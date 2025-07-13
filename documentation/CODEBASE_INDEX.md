# Common Dictionary Codebase Index

## Project Overview
The Common Dictionary project is a domain-specific entity extraction framework focused on chemistry data. It provides a modular architecture for extracting, processing, and managing chemical entities from various sources.

## Directory Structure

### Root Level
- `README.md` - Project documentation and setup instructions
- `requirements.txt` - Python dependencies
- `run_extractor_testing.sh` - Shell script for running extraction tests
- `.venv/` - Python virtual environment

### Core Directories

#### `/src` - Main Source Code
- `config_reconciliation.py` - Configuration management and reconciliation logic
- `/config/` - Configuration loading utilities
- `/domains/` - Domain-specific implementations

#### `/config` - Configuration Files
- `/domains/chemistry/` - Chemistry-specific configurations
  - `conflict_resolution.yaml` - Conflict resolution rules
  - `entity_config.yaml` - Entity extraction configuration
  - `extraction_config.yaml` - Extraction process configuration
  - `/subdomains/` - Subdomain-specific configurations
    - `/inorganic/` - Inorganic chemistry configurations
    - `/organic/` - Organic chemistry configurations

#### `/data` - Data Storage
- `/outputs/chemistry/` - Extraction outputs
  - `/subdomains/` - Subdomain-specific outputs
  - `entities.csv` - Extracted entities in CSV format
  - `entities.html` - Extracted entities in HTML format
  - `entities.json` - Extracted entities in JSON format

#### `/documentation` - Project Documentation
- Configuration design documents
- Architecture documentation
- API documentation

#### `/env_templates` - Environment Templates
- `/chemistry/` - Chemistry-specific environment templates
  - `/config_templates/` - Configuration templates
  - `env.development` - Development environment variables
  - `env.production` - Production environment variables
  - `env.staging` - Staging environment variables

#### `/interface` - User Interface
- `streamlit_app.py` - Streamlit web interface

## Key Source Files

### Core Framework Files

#### `src/config_reconciliation.py` (18KB, 400 lines)
- **Purpose**: Configuration management and reconciliation
- **Key Features**:
  - Loads and validates configuration files
  - Handles configuration conflicts
  - Manages environment-specific settings
  - Provides configuration validation

#### `src/config/config_loader.py`
- **Purpose**: Configuration file loading utilities
- **Key Features**:
  - YAML configuration parsing
  - Environment variable handling
  - Configuration validation

#### `src/config/env_loader.py`
- **Purpose**: Environment variable management
- **Key Features**:
  - Loads environment-specific configurations
  - Handles sensitive configuration data
  - Provides configuration overrides

### Chemistry Domain Files

#### `src/domains/chemistry/main.py` (2.3KB, 49 lines)
- **Purpose**: Main entry point for chemistry extraction
- **Key Features**:
  - Orchestrates extraction process
  - Manages data flow between components
  - Handles command-line interface

#### `src/domains/chemistry/extractor.py` (12KB, 225 lines)
- **Purpose**: Core extraction logic for chemistry entities
- **Key Features**:
  - Entity extraction from multiple sources
  - Data processing and validation
  - Result aggregation and deduplication

#### `src/domains/chemistry/output_generator.py` (11KB, 235 lines)
- **Purpose**: Output formatting and generation
- **Key Features**:
  - CSV, JSON, and HTML output generation
  - Data formatting and presentation
  - Export functionality

#### `src/domains/chemistry/example_usage.py` (6.3KB, 174 lines)
- **Purpose**: Usage examples and demonstrations
- **Key Features**:
  - Example extraction workflows
  - Configuration examples
  - Best practices demonstration

### Data Sources

#### `src/domains/chemistry/sources/base_source.py` (981B, 30 lines)
- **Purpose**: Abstract base class for data sources
- **Key Features**:
  - Defines source interface
  - Common source functionality
  - Error handling patterns

#### `src/domains/chemistry/sources/pubchem_source.py` (7.8KB, 194 lines)
- **Purpose**: PubChem API integration
- **Key Features**:
  - PubChem compound data extraction
  - Property mapping and transformation
  - Batch processing capabilities
  - Multiple output formats (JSON, CSV, TXT)

#### `src/domains/chemistry/sources/rdkit_source.py` (2.1KB, 62 lines)
- **Purpose**: RDKit molecular processing
- **Key Features**:
  - Molecular structure analysis
  - Chemical property calculations
  - Molecular descriptor generation

#### `src/domains/chemistry/sources/source_factory.py` (3.5KB, 91 lines)
- **Purpose**: Source instantiation and management
- **Key Features**:
  - Dynamic source creation
  - Source configuration management
  - Factory pattern implementation

#### `src/domains/chemistry/sources/test_source.py` (1.1KB, 34 lines)
- **Purpose**: Testing and validation source
- **Key Features**:
  - Mock data generation
  - Testing utilities
  - Validation helpers

## Configuration System

### Configuration Hierarchy
1. **Base Configuration**: Default settings in config files
2. **Environment Configuration**: Environment-specific overrides
3. **Domain Configuration**: Domain-specific settings
4. **Subdomain Configuration**: Subdomain-specific customizations

### Key Configuration Files
- `entity_config.yaml` - Entity extraction rules and mappings
- `extraction_config.yaml` - Extraction process parameters
- `conflict_resolution.yaml` - Data conflict resolution strategies

## Data Flow

### Extraction Pipeline
1. **Configuration Loading**: Load and validate configurations
2. **Source Initialization**: Initialize data sources
3. **Data Extraction**: Extract data from sources
4. **Processing**: Process and validate extracted data
5. **Conflict Resolution**: Resolve data conflicts
6. **Output Generation**: Generate formatted outputs

### Output Formats
- **CSV**: Tabular data for analysis
- **JSON**: Structured data for APIs
- **HTML**: Human-readable reports

## Environment Management

### Environment Types
- **Development**: Local development settings
- **Staging**: Pre-production testing
- **Production**: Live deployment settings

### Environment Variables
- Database connections
- API keys and credentials
- Logging configurations
- Performance settings

## Testing and Validation

### Testing Components
- Unit tests for individual components
- Integration tests for data flow
- End-to-end extraction tests
- Configuration validation tests

### Validation Features
- Data quality checks
- Configuration validation
- Output format validation
- Error handling and reporting

## Usage Examples

### Basic Extraction
```python
from src.domains.chemistry.main import run_extraction

# Run extraction with default configuration
results = run_extraction()
```

### Custom Configuration
```python
from src.domains.chemistry.extractor import ChemistryExtractor

# Initialize with custom configuration
extractor = ChemistryExtractor(config_path="custom_config.yaml")
results = extractor.extract()
```

### Source-Specific Extraction
```python
from src.domains.chemistry.sources.source_factory import SourceFactory

# Create specific source
pubchem_source = SourceFactory.create_source("pubchem")
data = pubchem_source.extract(cids=["172642621"])
```

## Dependencies

### Core Dependencies
- `requests` - HTTP client for API calls
- `pandas` - Data manipulation and analysis
- `pyyaml` - YAML configuration parsing
- `rdkit` - Chemical informatics toolkit

### Development Dependencies
- `pytest` - Testing framework
- `black` - Code formatting
- `flake8` - Code linting

## Architecture Patterns

### Design Patterns Used
1. **Factory Pattern**: Source creation and management
2. **Strategy Pattern**: Different extraction strategies
3. **Template Method**: Base source implementation
4. **Configuration Pattern**: Flexible configuration management

### Key Principles
- **Modularity**: Separate concerns into distinct modules
- **Extensibility**: Easy to add new sources and domains
- **Configurability**: Flexible configuration system
- **Testability**: Comprehensive testing support

## Performance Considerations

### Optimization Features
- Batch processing for API calls
- Caching mechanisms
- Parallel processing capabilities
- Memory-efficient data handling

### Monitoring and Logging
- Comprehensive logging system
- Performance metrics collection
- Error tracking and reporting
- Debug mode for troubleshooting

## Security Features

### Data Protection
- Secure credential management
- API key protection
- Data validation and sanitization
- Access control mechanisms

### Error Handling
- Graceful error recovery
- Detailed error reporting
- Fallback mechanisms
- Input validation

## Future Enhancements

### Planned Features
- Additional data sources
- Advanced conflict resolution
- Machine learning integration
- Real-time processing capabilities
- Enhanced UI/UX

### Scalability Improvements
- Distributed processing
- Database integration
- Caching layer
- API rate limiting 