# Technical Code Index - Common Dictionary

## Core Classes and Functions

### Configuration Management

#### `src/config_reconciliation.py`
```python
class ConfigReconciliation:
    def __init__(self, base_config_path: str, env_config_path: str = None)
    def load_configurations(self) -> Dict
    def reconcile_configs(self, base_config: Dict, env_config: Dict) -> Dict
    def validate_config(self, config: Dict) -> bool
    def merge_configs(self, configs: List[Dict]) -> Dict
```

#### `src/config/config_loader.py`
```python
class ConfigLoader:
    def __init__(self, config_dir: str)
    def load_yaml(self, file_path: str) -> Dict
    def load_config(self, config_name: str) -> Dict
    def validate_schema(self, config: Dict, schema: Dict) -> bool
```

#### `src/config/env_loader.py`
```python
class EnvLoader:
    def __init__(self, env_file: str)
    def load_env_vars(self) -> Dict[str, str]
    def get_config_value(self, key: str, default: Any = None) -> Any
    def set_env_var(self, key: str, value: str) -> None
```

### Chemistry Domain

#### `src/domains/chemistry/main.py`
```python
def run_extraction(config_path: str = None, output_format: str = "all") -> Dict
def setup_logging(log_level: str = "INFO") -> None
def validate_inputs(cids: List[str]) -> bool
```

#### `src/domains/chemistry/extractor.py`
```python
class ChemistryExtractor:
    def __init__(self, config: Dict)
    def extract(self, cids: List[str] = None) -> Dict
    def process_results(self, raw_data: Dict) -> Dict
    def resolve_conflicts(self, data: List[Dict]) -> List[Dict]
    def validate_entities(self, entities: List[Dict]) -> List[Dict]
```

#### `src/domains/chemistry/output_generator.py`
```python
class OutputGenerator:
    def __init__(self, output_dir: str)
    def generate_csv(self, data: List[Dict], filename: str) -> str
    def generate_json(self, data: List[Dict], filename: str) -> str
    def generate_html(self, data: List[Dict], filename: str) -> str
    def format_data(self, data: List[Dict], format_type: str) -> str
```

### Data Sources

#### `src/domains/chemistry/sources/base_source.py`
```python
class BaseSource(ABC):
    def __init__(self, config: Dict)
    @abstractmethod
    def extract(self, **kwargs) -> Dict
    def validate_response(self, response: Any) -> bool
    def handle_error(self, error: Exception) -> None
```

#### `src/domains/chemistry/sources/pubchem_source.py`
```python
class PubChemSource(BaseSource):
    def __init__(self, config: Dict)
    def extract(self, cids: List[str] = None) -> Dict
    def extract_properties(self, cids: List[str], properties: List[str]) -> Dict
    def extract_batch(self, cids: List[str], batch_size: int = 100) -> List[Dict]
    def format_property_url(self, cid: str, properties: List[str], format_type: str) -> str
    def parse_response(self, response: requests.Response, format_type: str) -> List[Dict]
```

#### `src/domains/chemistry/sources/rdkit_source.py`
```python
class RDKitSource(BaseSource):
    def __init__(self, config: Dict)
    def extract(self, smiles: List[str] = None) -> Dict
    def calculate_properties(self, mol: rdkit.Chem.Mol) -> Dict
    def generate_descriptors(self, mol: rdkit.Chem.Mol) -> Dict
```

#### `src/domains/chemistry/sources/source_factory.py`
```python
class SourceFactory:
    @staticmethod
    def create_source(source_type: str, config: Dict) -> BaseSource
    @staticmethod
    def get_available_sources() -> List[str]
    @staticmethod
    def validate_source_config(source_type: str, config: Dict) -> bool
```

## Configuration Schemas

### Entity Configuration (`entity_config.yaml`)
```yaml
entities:
  compound:
    properties:
      - molecular_weight
      - canonical_smiles
      - formula
      - iupac_name
    validation_rules:
      - required_fields: ["cid", "molecular_weight"]
      - data_types:
          cid: integer
          molecular_weight: float
          canonical_smiles: string
```

### Extraction Configuration (`extraction_config.yaml`)
```yaml
extraction:
  batch_size: 100
  max_retries: 3
  timeout: 30
  sources:
    pubchem:
      enabled: true
      rate_limit: 5
    rdkit:
      enabled: true
  output_formats:
    - csv
    - json
    - html
```

### Conflict Resolution (`conflict_resolution.yaml`)
```yaml
conflict_resolution:
  strategies:
    molecular_weight:
      priority: ["pubchem", "rdkit"]
      aggregation: "mean"
    canonical_smiles:
      priority: ["pubchem", "rdkit"]
      aggregation: "most_frequent"
```

## Data Structures

### Entity Data Structure
```python
Entity = {
    "cid": int,                    # PubChem Compound ID
    "molecular_weight": float,     # Molecular weight in g/mol
    "canonical_smiles": str,       # Canonical SMILES string
    "formula": str,               # Molecular formula
    "iupac_name": str,            # IUPAC name
    "source": str,                # Data source (pubchem, rdkit)
    "confidence": float,          # Confidence score (0-1)
    "timestamp": datetime,        # Extraction timestamp
    "metadata": Dict              # Additional metadata
}
```

### Extraction Result Structure
```python
ExtractionResult = {
    "entities": List[Entity],     # List of extracted entities
    "statistics": {
        "total_entities": int,
        "successful_extractions": int,
        "failed_extractions": int,
        "processing_time": float
    },
    "errors": List[str],          # List of error messages
    "warnings": List[str],        # List of warning messages
    "metadata": Dict              # Additional metadata
}
```

## API Endpoints and Data Sources

### PubChem API Endpoints
```python
# Property extraction
"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/{properties}/{format}"

# Batch property extraction
"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cids}/property/{properties}/{format}"

# Compound information
"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/JSON"
```

### Supported Properties
- `MolecularWeight` - Molecular weight in g/mol
- `CanonicalSMILES` - Canonical SMILES representation
- `Formula` - Molecular formula
- `IUPACName` - IUPAC systematic name
- `InChI` - International Chemical Identifier
- `InChIKey` - InChI key

### Output Formats
- `JSON` - JSON format (default)
- `CSV` - Comma-separated values
- `TXT` - Plain text format

## Error Handling Patterns

### Source Error Handling
```python
def handle_source_error(self, error: Exception, context: Dict) -> None:
    """Handle errors from data sources"""
    error_type = type(error).__name__
    error_msg = str(error)
    
    if isinstance(error, requests.RequestException):
        self.logger.error(f"Network error: {error_msg}")
        # Implement retry logic
    elif isinstance(error, ValueError):
        self.logger.error(f"Data validation error: {error_msg}")
        # Handle invalid data
    else:
        self.logger.error(f"Unexpected error: {error_msg}")
        # Handle unknown errors
```

### Configuration Error Handling
```python
def validate_config(self, config: Dict) -> bool:
    """Validate configuration structure and values"""
    required_fields = ["entities", "extraction", "sources"]
    
    for field in required_fields:
        if field not in config:
            raise ValueError(f"Missing required field: {field}")
    
    # Validate nested structures
    self._validate_entities_config(config["entities"])
    self._validate_extraction_config(config["extraction"])
    
    return True
```

## Logging Configuration

### Logging Setup
```python
def setup_logging(log_level: str = "INFO", log_file: str = None) -> None:
    """Setup logging configuration"""
    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(log_file) if log_file else logging.NullHandler()
        ]
    )
```

### Debug Logging
```python
def debug_log(self, message: str, data: Any = None) -> None:
    """Log debug information with optional data"""
    if self.debug_enabled:
        self.logger.debug(message)
        if data:
            self.logger.debug(f"Data: {data}")
```

## Performance Optimizations

### Batch Processing
```python
def process_batch(self, items: List[Any], batch_size: int) -> List[Any]:
    """Process items in batches for better performance"""
    results = []
    for i in range(0, len(items), batch_size):
        batch = items[i:i + batch_size]
        batch_results = self._process_single_batch(batch)
        results.extend(batch_results)
    return results
```

### Caching Implementation
```python
class CacheManager:
    def __init__(self, cache_dir: str, max_size: int = 1000):
        self.cache_dir = cache_dir
        self.max_size = max_size
        self.cache = {}
    
    def get(self, key: str) -> Any:
        """Get cached value"""
        if key in self.cache:
            return self.cache[key]
        return None
    
    def set(self, key: str, value: Any) -> None:
        """Set cached value"""
        if len(self.cache) >= self.max_size:
            self._evict_oldest()
        self.cache[key] = value
```

## Testing Patterns

### Unit Test Structure
```python
class TestPubChemSource(unittest.TestCase):
    def setUp(self):
        """Setup test fixtures"""
        self.config = {
            "base_url": "https://pubchem.ncbi.nlm.nih.gov/rest/pug",
            "timeout": 30,
            "max_retries": 3
        }
        self.source = PubChemSource(self.config)
    
    def test_extract_single_cid(self):
        """Test extraction of single CID"""
        result = self.source.extract(cids=["172642621"])
        self.assertIsInstance(result, dict)
        self.assertIn("entities", result)
    
    def test_extract_batch(self):
        """Test batch extraction"""
        cids = ["172642621", "172642622"]
        result = self.source.extract_batch(cids, batch_size=2)
        self.assertEqual(len(result), 2)
```

### Integration Test Structure
```python
class TestChemistryExtractor(unittest.TestCase):
    def setUp(self):
        """Setup integration test fixtures"""
        self.config = self.load_test_config()
        self.extractor = ChemistryExtractor(self.config)
    
    def test_full_extraction_pipeline(self):
        """Test complete extraction pipeline"""
        cids = ["172642621"]
        result = self.extractor.extract(cids)
        
        # Validate result structure
        self.assertIn("entities", result)
        self.assertIn("statistics", result)
        self.assertGreater(len(result["entities"]), 0)
```

## Environment Variables

### Required Environment Variables
```bash
# Database Configuration
DATABASE_URL=postgresql://user:pass@localhost/dbname

# API Configuration
PUBCHEM_API_KEY=your_api_key_here
PUBCHEM_RATE_LIMIT=5

# Logging Configuration
LOG_LEVEL=INFO
LOG_FILE=extraction.log

# Performance Configuration
BATCH_SIZE=100
MAX_RETRIES=3
TIMEOUT=30

# Debug Configuration
DEBUG_PUBCHEM=true
DEBUG_RDKIT=false
```

### Optional Environment Variables
```bash
# Output Configuration
OUTPUT_DIR=/path/to/output
OUTPUT_FORMATS=csv,json,html

# Cache Configuration
CACHE_ENABLED=true
CACHE_DIR=/path/to/cache
CACHE_MAX_SIZE=1000

# Monitoring Configuration
METRICS_ENABLED=true
METRICS_PORT=8080
```

## Deployment Configuration

### Docker Configuration
```dockerfile
FROM python:3.9-slim

WORKDIR /app

COPY requirements.txt .
RUN pip install -r requirements.txt

COPY . .

CMD ["python", "-m", "src.domains.chemistry.main"]
```

### Docker Compose
```yaml
version: '3.8'
services:
  extractor:
    build: .
    environment:
      - DATABASE_URL=postgresql://user:pass@db:5432/extractor
      - LOG_LEVEL=INFO
    volumes:
      - ./data:/app/data
      - ./logs:/app/logs
    depends_on:
      - db
  
  db:
    image: postgres:13
    environment:
      - POSTGRES_DB=extractor
      - POSTGRES_USER=user
      - POSTGRES_PASSWORD=pass
    volumes:
      - postgres_data:/var/lib/postgresql/data

volumes:
  postgres_data:
```

## Monitoring and Metrics

### Performance Metrics
```python
class MetricsCollector:
    def __init__(self):
        self.metrics = {
            "extraction_count": 0,
            "success_count": 0,
            "error_count": 0,
            "processing_time": 0.0
        }
    
    def record_extraction(self, success: bool, processing_time: float):
        """Record extraction metrics"""
        self.metrics["extraction_count"] += 1
        if success:
            self.metrics["success_count"] += 1
        else:
            self.metrics["error_count"] += 1
        self.metrics["processing_time"] += processing_time
    
    def get_metrics(self) -> Dict:
        """Get current metrics"""
        return self.metrics.copy()
```

### Health Check Endpoint
```python
def health_check() -> Dict:
    """Health check endpoint for monitoring"""
    return {
        "status": "healthy",
        "timestamp": datetime.now().isoformat(),
        "version": "1.0.0",
        "uptime": get_uptime(),
        "memory_usage": get_memory_usage(),
        "disk_usage": get_disk_usage()
    }
``` 