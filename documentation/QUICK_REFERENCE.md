# Quick Reference Guide - Common Dictionary

## Common Operations

### Running the Extractor

#### Basic Extraction
```bash
# Run with default configuration
cd common_dictionary
python -m src.domains.chemistry.main

# Run with custom configuration
python -m src.domains.chemistry.main --config config/domains/chemistry/extraction_config.yaml

# Run with specific CIDs
python -m src.domains.chemistry.main --cids 172642621,172642622
```

#### Using the Shell Script
```bash
# Run extraction testing
./run_extractor_testing.sh

# Run with debug mode
DEBUG_PUBCHEM=true ./run_extractor_testing.sh
```

### Configuration Management

#### Environment Setup
```bash
# Copy environment template
cp env_templates/chemistry/env.development .env

# Edit environment variables
nano .env

# Load environment
source .env
```

#### Configuration Validation
```bash
# Validate configuration files
python -c "
from src.config.config_loader import ConfigLoader
loader = ConfigLoader('config')
config = loader.load_config('domains/chemistry/entity_config.yaml')
print('Configuration is valid')
"
```

### Data Source Operations

#### PubChem Source
```python
from src.domains.chemistry.sources.pubchem_source import PubChemSource

# Initialize source
config = {"base_url": "https://pubchem.ncbi.nlm.nih.gov/rest/pug"}
source = PubChemSource(config)

# Extract single CID
result = source.extract(cids=["172642621"])

# Extract multiple properties
properties = ["MolecularWeight", "CanonicalSMILES", "Formula"]
result = source.extract_properties(["172642621"], properties)
```

#### RDKit Source
```python
from src.domains.chemistry.sources.rdkit_source import RDKitSource

# Initialize source
config = {}
source = RDKitSource(config)

# Extract from SMILES
smiles_list = ["CCO", "CC(=O)O"]
result = source.extract(smiles=smiles_list)
```

### Output Generation

#### Generate Different Formats
```python
from src.domains.chemistry.output_generator import OutputGenerator

# Initialize generator
generator = OutputGenerator("data/outputs/chemistry")

# Generate outputs
data = [{"cid": 172642621, "molecular_weight": 180.16}]

# CSV output
generator.generate_csv(data, "entities.csv")

# JSON output
generator.generate_json(data, "entities.json")

# HTML output
generator.generate_html(data, "entities.html")
```

## Troubleshooting

### Common Issues

#### 1. PubChem API Errors
**Problem**: Getting 404 or 500 errors from PubChem API
```bash
# Check API endpoint manually
curl "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/172642621/property/MolecularWeight/JSON"

# Verify CID exists
curl "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/172642621/JSON"
```

**Solutions**:
- Verify CID exists in PubChem
- Check property names are correct
- Ensure proper URL formatting
- Add rate limiting if hitting API limits

#### 2. Configuration Errors
**Problem**: Configuration validation fails
```bash
# Check configuration syntax
python -c "
import yaml
with open('config/domains/chemistry/entity_config.yaml') as f:
    yaml.safe_load(f)
print('YAML syntax is valid')
"
```

**Solutions**:
- Validate YAML syntax
- Check required fields are present
- Verify file paths are correct
- Ensure environment variables are set

#### 3. Empty Results
**Problem**: Extraction returns empty results
```python
# Enable debug logging
import logging
logging.basicConfig(level=logging.DEBUG)

# Check source responses
source = PubChemSource(config)
response = source.extract(cids=["172642621"])
print(f"Response: {response}")
```

**Solutions**:
- Enable debug mode (`DEBUG_PUBCHEM=true`)
- Check API responses manually
- Verify CIDs are valid
- Check property mappings

#### 4. Memory Issues
**Problem**: High memory usage with large datasets
```python
# Use batch processing
extractor = ChemistryExtractor(config)
result = extractor.extract_batch(cids, batch_size=50)
```

**Solutions**:
- Reduce batch size
- Process data in chunks
- Use streaming for large files
- Monitor memory usage

### Debug Commands

#### Enable Debug Logging
```bash
# Set debug environment variables
export DEBUG_PUBCHEM=true
export DEBUG_RDKIT=true
export LOG_LEVEL=DEBUG

# Run with debug output
python -m src.domains.chemistry.main
```

#### Check Logs
```bash
# View recent logs
tail -f data/outputs/chemistry/logs/chemistry_extraction.log

# Search for errors
grep -i error data/outputs/chemistry/logs/chemistry_extraction.log

# Search for specific CID
grep "172642621" data/outputs/chemistry/logs/chemistry_extraction.log
```

#### Validate Data
```python
# Check output files
import pandas as pd
import json

# CSV validation
df = pd.read_csv("data/outputs/chemistry/entities.csv")
print(f"CSV rows: {len(df)}")
print(f"CSV columns: {df.columns.tolist()}")

# JSON validation
with open("data/outputs/chemistry/entities.json") as f:
    data = json.load(f)
print(f"JSON entities: {len(data.get('entities', []))}")
```

### Performance Optimization

#### Batch Processing
```python
# Optimize batch size
config = {
    "extraction": {
        "batch_size": 100,  # Adjust based on API limits
        "max_retries": 3,
        "timeout": 30
    }
}
```

#### Caching
```python
# Enable caching
config = {
    "cache": {
        "enabled": True,
        "cache_dir": "cache",
        "max_size": 1000
    }
}
```

#### Rate Limiting
```python
# Configure rate limiting
config = {
    "sources": {
        "pubchem": {
            "rate_limit": 5,  # requests per second
            "delay": 0.2      # seconds between requests
        }
    }
}
```

## API Reference

### PubChem API Properties
| Property | Format | Description |
|----------|--------|-------------|
| `MolecularWeight` | JSON/CSV/TXT | Molecular weight in g/mol |
| `CanonicalSMILES` | JSON/CSV/TXT | Canonical SMILES string |
| `Formula` | JSON/CSV/TXT | Molecular formula |
| `IUPACName` | JSON/CSV/TXT | IUPAC systematic name |
| `InChI` | JSON/CSV/TXT | International Chemical Identifier |
| `InChIKey` | JSON/CSV/TXT | InChI key |

### Configuration Keys
| Key | Type | Description | Default |
|-----|------|-------------|---------|
| `batch_size` | int | Number of items per batch | 100 |
| `max_retries` | int | Maximum retry attempts | 3 |
| `timeout` | int | Request timeout in seconds | 30 |
| `rate_limit` | int | Requests per second | 5 |
| `debug_enabled` | bool | Enable debug logging | false |

### Environment Variables
| Variable | Description | Example |
|----------|-------------|---------|
| `DEBUG_PUBCHEM` | Enable PubChem debug logging | `true` |
| `DEBUG_RDKIT` | Enable RDKit debug logging | `false` |
| `LOG_LEVEL` | Logging level | `INFO` |
| `BATCH_SIZE` | Default batch size | `100` |
| `PUBCHEM_RATE_LIMIT` | PubChem API rate limit | `5` |

## File Locations

### Configuration Files
- `config/domains/chemistry/entity_config.yaml` - Entity extraction rules
- `config/domains/chemistry/extraction_config.yaml` - Extraction parameters
- `config/domains/chemistry/conflict_resolution.yaml` - Conflict resolution

### Output Files
- `data/outputs/chemistry/entities.csv` - CSV output
- `data/outputs/chemistry/entities.json` - JSON output
- `data/outputs/chemistry/entities.html` - HTML output
- `data/outputs/chemistry/logs/` - Log files

### Environment Files
- `env_templates/chemistry/env.development` - Development environment
- `env_templates/chemistry/env.production` - Production environment
- `env_templates/chemistry/env.staging` - Staging environment

## Testing

### Unit Tests
```bash
# Run all tests
python -m pytest tests/

# Run specific test file
python -m pytest tests/test_pubchem_source.py

# Run with verbose output
python -m pytest -v tests/
```

### Integration Tests
```bash
# Run integration tests
python -m pytest tests/integration/

# Run with coverage
python -m pytest --cov=src tests/
```

### Manual Testing
```python
# Test PubChem source
from src.domains.chemistry.sources.pubchem_source import PubChemSource
source = PubChemSource({})
result = source.extract(cids=["172642621"])
print(result)

# Test RDKit source
from src.domains.chemistry.sources.rdkit_source import RDKitSource
source = RDKitSource({})
result = source.extract(smiles=["CCO"])
print(result)
```

## Monitoring

### Health Checks
```python
# Check system health
def health_check():
    return {
        "status": "healthy",
        "timestamp": datetime.now().isoformat(),
        "memory_usage": psutil.virtual_memory().percent,
        "disk_usage": psutil.disk_usage('/').percent
    }
```

### Metrics Collection
```python
# Collect performance metrics
from src.utils.metrics import MetricsCollector
metrics = MetricsCollector()
metrics.record_extraction(success=True, processing_time=1.5)
print(metrics.get_metrics())
```

## Deployment

### Docker Deployment
```bash
# Build image
docker build -t common-dictionary .

# Run container
docker run -d \
  -e DEBUG_PUBCHEM=true \
  -v $(pwd)/data:/app/data \
  common-dictionary

# Run with docker-compose
docker-compose up -d
```

### Production Deployment
```bash
# Set production environment
cp env_templates/chemistry/env.production .env

# Install dependencies
pip install -r requirements.txt

# Run with production settings
python -m src.domains.chemistry.main
``` 