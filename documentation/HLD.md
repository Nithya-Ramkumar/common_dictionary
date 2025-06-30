# Common Dictionary Module — Chemistry Domain: Executive Summary

The Common Dictionary Module is a modular, extensible, and provenance-aware system for managing chemical knowledge. It is designed to support robust extraction, enrichment, validation, and human-in-the-loop review of chemical entities and relationships, with a focus on traceability, data quality, and industry best practices.

**Key Design Principles:**
- The core configuration (`entity_config.yaml`) defines only the schema (structure, attributes, validation rules) for each entity type (e.g., Compound, Reaction, Catalyst, Solvent). It does not contain any real-world instances or static UUIDs.
- UUIDs are not present in the schema config. Instead, they are generated and assigned to each unique, reconciled entity instance only after extraction, reconciliation, and generation. This ensures global uniqueness, traceability, and avoids manual errors.
- Clear separation of concerns: schema definition, extraction configuration, source priority mapping, conflict resolution, and validation are all handled in separate, modular YAML files.
- **Environment-aware configuration system** with support for multiple environments (development, production, testing, staging) and automatic environment detection.
- For every field, all source values, alternatives, and their provenance are tracked—not just the canonical value. This supports transparency, auditability, and robust human review.
- A Streamlit UI presents all source values, alternatives, provenance, and occurrences for each field. Reviewers can approve, override, or edit values, and all actions are logged. The reviewed and reconciled output becomes the single source of truth for downstream use.
- All quantitative fields are paired with value and unit fields. Units are validated against a central, shared `metric_units.yaml` (generated from the Pint library and SI/NIST standards), ensuring consistency and extensibility.
- A reconciliation script ensures that every entity field has a source mapping and conflict resolution strategy, and that there are no orphaned mappings. The system will not proceed unless all configs are consistent.

**Workflow Summary:**
1. **Environment setup** - Load environment-specific configuration (API keys, database settings, feature flags, etc.)
2. Define entity types and their attributes in `entity_config.yaml` (no UUIDs, no instances).
3. Extract raw data from multiple sources using the extraction config and environment-specific settings.
4. Merge, deduplicate, and enrich extracted data. Assign UUIDs to each unique, reconciled entity instance at this stage.
5. Validate all data against the schema and validation config, including value ranges, regex, cross-field logic, and synonym clustering.
6. Apply conflict resolution strategies as defined in the config, logging all alternatives and provenance.
7. Present all data, provenance, and alternatives in a UI for human review and approval. Log all reviewer actions.
8. Produce a `.json` output containing all reviewed, validated, and UUID-assigned entity instances. This output is the single source of truth for downstream extraction and analysis.

**Output and Traceability:**
- The final output contains, for each entity instance: a UUID (assigned post-extraction), all attributes and their selected values, all alternative values and their provenance, full review and edit history, and versioning/timestamps.
- This output ensures full traceability, auditability, and compliance with industry standards.

**Extensibility:**
- The design supports easy addition of new domains, subdomains, entity types, and attributes. All configs are modular and can be overridden or extended at the subdomain level.
- **Environment configuration** supports easy addition of new environments and environment-specific settings.

**Summary Table:**
| Aspect                | Approach/Best Practice                                      |
|-----------------------|------------------------------------------------------------|
| UUIDs                 | Assigned only to instances post-extraction, not in schema  |
| Schema Config         | Defines structure, attributes, and rules (no instances)    |
| Extraction Config     | Source connection details and parameters only              |
| Source Priority       | Attribute-to-source mapping with ordered priorities        |
| Conflict Resolution   | Config-driven, entity-wise, and relationship-wise          |
| Validation            | Schema, value, regex, cross-field, and synonym validation  |
| Units                 | Centralized, SI-compliant, validated via shared config     |
| Provenance            | All source values and alternatives tracked                 |
| Human Review          | UI for review, approval, and audit trail                   |
| Output                | JSON with UUIDs, provenance, review history, versioning    |
| **Environment Config**| **Multi-environment support with automatic detection**     |
vscode
**In summary:**
The Common Dictionary Module is a robust, industry-aligned system for managing chemical knowledge, ensuring modularity, traceability, provenance, and human oversight at every stage. Its design is future-proof, extensible, and ready for integration with downstream scientific workflows.

---

# Common Dictionary Module: Chemistry Domain — Design Document

---

## 1. Purpose

To create a **modular, extensible, and provenance-aware dictionary system** for the chemistry domain, supporting:
- Consistent, schema-validated entity/relationship/graph definitions
- Multi-source, priority-driven enrichment with full provenance tracking
- UUID-based traceability for all objects (assigned only to instances post-extraction)
- **Environment-aware configuration** with support for multiple deployment environments
- Human-in-the-loop review and iterative improvement
- Exportable, versioned `.json` output for downstream extraction from publications and journals

---

## 2. Key Design Principles

This section expands on the principles outlined in the executive summary:
- **UUIDs:** Only assigned to entity instances after extraction and reconciliation, never in the schema config.
- **Separation of Concerns:** Each config file has a single responsibility (schema, extraction, source priority, conflict resolution, validation).
- **Provenance:** All source values and alternatives are tracked for every field.
- **Human-in-the-Loop:** UI-driven review and approval, with full audit trail.
- **Extensibility:** Modular configs and code structure for easy domain and subdomain expansion.
- **Industry Standards:** SI-compliant units, schema validation, and best practices throughout.
- **Environment Management:** Multi-environment support with automatic detection and configuration loading.

---

## 3. Environment Configuration System

### 3.1 Overview

The system includes a comprehensive environment configuration system that supports multiple deployment environments (development, production, testing, staging) with automatic environment detection and configuration loading.

### 3.2 Environment Structure

```
common_dictionary/
├── env_templates/
│   └── chemistry/
│       ├── env.template          # Comprehensive template with all variables
│       ├── env.development       # Development-specific settings
│       ├── env.production        # Production-specific settings
│       ├── env.testing           # Testing-specific settings
│       └── env.staging           # Staging-specific settings
```

### 3.3 Environment Variables

The environment system manages the following categories of configuration:

#### **Chemistry Domain Configuration Paths**
- `CONFIG_ROOT` - Root configuration directory
- `CONFIG_TEST` - Test configuration directory
- Individual config file paths (entity, extraction, source mapping, etc.)

#### **Chemistry Data Source API Settings**
- **PubChem API:** Base URL, timeout, retry settings, batch size
- **Reaxys API:** Base URL, API key, timeout, retry settings
- **ChEBI API:** Base URL, timeout, retry settings

#### **Rate Limiting & Retry Configuration**
- `MAX_REQUESTS_PER_MINUTE` - API rate limiting
- `RATE_LIMIT_WINDOW` - Rate limit window in seconds
- `RETRY_ATTEMPTS` - Number of retry attempts
- `RETRY_DELAY` - Delay between retries
- `MAX_CONCURRENT_DOWNLOADS` - Concurrent download limit

#### **Storage Configuration**
- `DATA_ROOT_DIR` - Root data directory
- Chemistry-specific storage paths (raw, processed, metrics)
- Document and output directories

#### **Logging Configuration**
- `LOG_LEVEL` - Logging level (DEBUG, INFO, WARNING, ERROR)
- `LOG_FILE` - Log file path
- `LOG_FORMAT` - Log message format
- `LOG_ROTATION` - Log rotation settings
- `LOG_RETENTION` - Log retention period

#### **Database Configuration**
- `DB_HOST`, `DB_PORT`, `DB_NAME` - Database connection settings
- `DB_USER`, `DB_PASSWORD` - Database credentials
- `DB_POOL_SIZE`, `DB_MAX_OVERFLOW` - Connection pool settings

#### **Cache Configuration**
- `CACHE_TYPE` - Cache type (redis, memory)
- `CACHE_HOST`, `CACHE_PORT` - Cache server settings
- `CACHE_TTL` - Cache time-to-live
- `CACHE_MAX_SIZE` - Maximum cache size
- `CACHE_ENABLED` - Cache enable/disable flag

#### **Security Settings**
- `API_KEY_EXPIRY_DAYS` - API key expiration
- `JWT_SECRET_KEY` - JWT secret key
- `JWT_ALGORITHM` - JWT algorithm
- `TOKEN_EXPIRY_MINUTES` - Token expiration time
- `HASH_ALGORITHM` - Password hashing algorithm

#### **Proxy Configuration**
- `USE_PROXY` - Enable/disable proxy
- `HTTP_PROXY`, `HTTPS_PROXY` - Proxy server settings
- `NO_PROXY` - Proxy exclusion list

#### **Error Handling**
- `DETAILED_ERRORS` - Enable detailed error messages
- `ALERT_ON_ERROR` - Enable error alerts
- `ALERT_EMAIL` - Alert email address

#### **Feature Flags**
- `ENABLE_RATE_LIMITING` - Enable/disable rate limiting
- `ENABLE_CACHING` - Enable/disable caching
- `ENABLE_MONITORING` - Enable/disable monitoring
- `ENABLE_ANALYTICS` - Enable/disable analytics

#### **Environment Settings**
- `ENVIRONMENT` - Current environment (development, production, testing, staging)
- `DEBUG` - Debug mode flag

#### **Chemistry-Specific Features**
- `ENABLE_PUBCHEM`, `ENABLE_REAXYS`, `ENABLE_CHEBI` - Enable/disable data sources
- `CHEMISTRY_BATCH_PROCESSING` - Enable batch processing
- `CHEMISTRY_VALIDATION_STRICT` - Enable strict validation
- `CHEMISTRY_SYNONYM_CLUSTERING` - Enable synonym clustering
- Default units for mass, volume, temperature, pressure

### 3.4 Environment Loader

The `EnvironmentLoader` class provides:

- **Automatic environment detection** based on `ENVIRONMENT` variable
- **Fallback mechanism** to template if environment-specific file not found
- **Type conversion** and validation for all environment variables
- **Specialized getter methods** for different configuration sections
- **Environment status checking** (is_development, is_production, etc.)

### 3.5 Environment-Specific Optimizations

#### **Development Environment**
- Debug mode enabled
- Detailed error messages
- Lower rate limits and batch sizes
- Memory-based caching
- Extended token expiry for testing

#### **Production Environment**
- High performance settings
- Strict validation enabled
- All data sources enabled
- Redis caching
- Standard security settings

#### **Testing Environment**
- Minimal resource usage
- Mock data sources
- Disabled monitoring and analytics
- Short timeouts and retry limits
- Memory-based caching

#### **Staging Environment**
- Moderate settings for testing
- Limited data sources
- Debugging enabled
- Intermediate rate limits
- Shorter cache TTL

---

## 4. Entity & Relationship Definition

### 4.1 Entity Example: Compound

```yaml
- name: Compound
  description: "Chemical compounds and molecules"
  domain: chemistry
  attributes:
    - name: name
      type: string
      required: true
      validation:
        min_length: 2
        max_length: 100
      provenance: {}
    - name: formula
      type: string
      required: true
      validation:
        pattern: "^[A-Z][a-z]?[0-9]*$"
      provenance: {}
    # ... more attributes ...
  validation_rules:
    - rule: "required_fields"
      fields: ["name", "formula"]
    - rule: "formula_format"
      pattern: "^[A-Z][a-z]?[0-9]*$"
  enrichment_rules:
    - name: "pubchem_lookup"
      description: "Look up compound information from PubChem"
      source: "pubchem_api"
      enabled: true
    - name: "structure_parsing"
      description: "Parse SMILES structure"
      source: "rdkit"
      enabled: true
  review_history: []
  version: "1.0.0"
  last_modified: "2024-06-01T12:00:00Z"
```

> **Note:** No UUID is present in the schema/class config. UUIDs are assigned to specific compound instances only after extraction, reconciliation, and generation, and are present in the output files.

### 4.2 Relationship Example

```yaml
- name: catalyzes
  description: "Catalyst catalyzes reaction"
  domain: chemistry
  source_priority:
    - source: reaxys
      priority: 1
    - source: pubchem
      priority: 2
  provenance:
    selected_value: "Catalyst X catalyzes Reaction Y"
    selected_source: "reaxys"
    alternatives:
      - value: "Catalyst X catalyzes Reaction Y"
        source: "pubchem"
        confidence: 0.85
  validation_rules:
    - rule: "valid_types"
      enabled: true
  review_history: []
  version: "1.0.0"
  last_modified: "2024-06-01T12:00:00Z"
```

---

## 5. Source Priority Mapping

- The `source_priority.yaml` file contains, for each entity type and attribute, the ordered list of sources (with priority).
- Every entity field must have a priority and source mapping here.
- This mapping is used to select the canonical value and track all alternatives for provenance and human review.

---

## 6. Extraction Configuration

- The `extraction_config.yaml` file only contains source connection details, API endpoints, and general extraction parameters.
- It does not specify which attributes to extract or their priorities.
- Includes a `record_all_occurrences` flag to control whether all found locations/contexts are tracked for provenance.
- This separation allows for independent management of extraction logic and source prioritization.

---

## 7. Conflict Resolution

- The `conflict_resolution.yaml` file defines how to resolve conflicts when multiple sources provide different values for the same attribute.
- Entity-wise conflict resolution is currently populated (e.g., for Compound, Reaction, etc.).
- Relationship-wise conflict resolution is reserved for future expansion.
- Strategies can include: highest-priority source, highest confidence, human review, or custom logic.

---

## 8. Configuration Reconciliation

- A Python script (`config_reconciliation.py`) validates that:
  - Every attribute defined in `entity_config.yaml` for every entity type has a corresponding entry in `source_priority.yaml`.
  - There are no orphaned mappings in `source_priority.yaml` that do not correspond to any entity attribute.
  - All required conflict resolution strategies are defined for each entity.
- This validation is required before the system can proceed. If validation fails, the system will not run until configs are consistent.
- This ensures no attribute is left without a provenance/source mapping and keeps the system robust and auditable.

---

## 9. Enhanced Configuration Management

### 9.1 EnvironmentLoader Integration

The system now includes an `EnvironmentLoader` class that:

- **Automatically detects** the current environment
- **Loads environment-specific** configuration files
- **Provides type-safe access** to all environment variables
- **Supports fallback** to template configuration
- **Integrates seamlessly** with existing configuration loaders

### 9.2 Enhanced ConfigLoader

The `ConfigLoader` class has been enhanced to:

- **Integrate with EnvironmentLoader** for seamless environment management
- **Maintain backward compatibility** with existing code
- **Provide enhanced logging** based on environment configuration
- **Support comprehensive configuration** loading with error handling
- **Include environment-specific methods** for checking current environment
- **Provide access to storage paths** and rate limiting configuration

### 9.3 Source Factory Enhancement

The `SourceFactory` has been updated to:

- **Use EnvironmentLoader** for configuration
- **Support environment-aware** source creation
- **Validate source availability** based on environment settings
- **Provide enhanced error handling** and logging
- **Support multiple data sources** with environment-specific enablement

---

## 10. Extraction, Enrichment, and Generation Models/Logic

| Task                  | Model/Tool/Logic         | Notes/Source Priority         |
|-----------------------|-------------------------|------------------------------|
| NER Extraction        | SciSpacy, BERN, BioBERT | Domain-specific, configurable|
| Ontology Linking      | BioSyn, MetaMap         | ChEBI, PubChem, Reaxys       |
| Synonym Clustering    | SBERT, DBSCAN           | Semantic similarity          |
| SMILES Parsing        | RDKit, DeepChem         | Structure validation         |
| Enrichment            | PubChem, Reaxys, ChEBI  | API, priority order          |
| Graph Construction    | NetworkX                | Source-priority for edges    |
| Validation            | JSON Schema, Regex      | Config-driven                |
| UUID Assignment       | uuid4                   | For all instances, post-extraction |
| Review/Export         | Streamlit, Pandas       | Human-in-the-loop, exportable|
| **Environment Config**| **EnvironmentLoader**   | **Multi-environment support**|

---

## 11. Human-in-the-Loop Workflow

- Streamlit UI displays for each field:
  - The selected value (from the highest-priority source or as resolved by conflict strategy)
  - All alternative values, sources, and confidence
  - All occurrences/contexts where the value was found (if `record_all_occurrences` is enabled)
  - Reviewer can approve, override, or edit
  - All actions are logged (reviewer, timestamp, rationale)
- Iterative review: Each review/approval updates the `.json` output.
- Export: The `.json` output (with full provenance, all occurrences, review history, and UUIDs assigned to instances) is used for downstream extraction from publications and journals.

---

## 12. Output Format Example

```json
{
  "uuid": "d58e6e10-89c2-4d6e-bad0-9ecb9be9a8a7",
  "name": {
    "selected_value": "Aspirin",
    "selected_source": "pubchem",
    "alternatives": [
      {"value": "Acetylsalicylic acid", "source": "chebi", "confidence": 0.95},
      {"value": "2-Acetoxybenzoic acid", "source": "reaxys", "confidence": 0.92}
    ],
    "occurrences": [
      {"source": "pubchem", "context": "CID:2244"},
      {"source": "chebi", "context": "CHEBI:15365"},
      {"source": "reaxys", "context": "ReaxysID:12345"}
    ]
  },
  "formula": {
    "selected_value": "C9H8O4",
    "selected_source": "pubchem",
    "alternatives": [
      {"value": "C9H8O4", "source": "chebi", "confidence": 0.99}
    ]
  },
  "domain": "chemistry",
  "attributes": { ... },
  "validation_rules": [ ... ],
  "enrichment_rules": [ ... ],
  "review_history": [
    {
      "reviewer": "user1",
      "timestamp": "2024-06-01T12:00:00Z",
      "action": "approved",
      "field": "name",
      "selected_value": "Aspirin",
      "rationale": "PubChem is most authoritative"
    }
  ],
  "version": "1.0.0",
  "last_modified": "2024-06-01T12:00:00Z"
}
```

> **Note:** UUIDs are present only in the output files for specific instances, not in the schema config.

---

## 13. Updated Folder Structure

```
common_dictionary/
├── env_templates/                    # NEW: Environment configuration templates
│   └── chemistry/
│       ├── env.template              # Comprehensive template with all variables
│       ├── env.development           # Development-specific settings
│       ├── env.production            # Production-specific settings
│       ├── env.testing               # Testing-specific settings
│       └── env.staging               # Staging-specific settings
├── config/
│   └── domains/
│       └── chemistry/
│           ├── entity_config.yaml
│           ├── extraction_config.yaml
│           ├── validation_config.yaml
│           ├── source_priority.yaml
│           ├── conflict_resolution.yaml
│           ├── relationship_config.yaml
│           └── subdomains/
│               ├── organic/
│               │   ├── entity_config.yaml
│               │   ├── extraction_config.yaml
│               │   └── validation_config.yaml
│               ├── inorganic/
│               │   ├── entity_config.yaml
│               │   ├── extraction_config.yaml
│               │   └── validation_config.yaml
│               └── ... (other subdomains)
├── src/
│   ├── config/                       # NEW: Enhanced configuration management
│   │   ├── env_loader.py             # EnvironmentLoader class
│   │   └── config_loader.py          # Enhanced ConfigLoader class
│   └── domains/
│       └── chemistry/
│           ├── sources/              # Enhanced source management
│           │   ├── base_source.py
│           │   ├── source_factory.py # Updated with environment support
│           │   ├── pubchem_source.py # Updated with environment support
│           │   └── test_source.py
│           ├── extractors/
│           ├── processors/
│           ├── generators/
│           ├── utils/
│           ├── config_reconciliation.py
│           ├── example_usage.py      # Updated with environment examples
│           └── subdomains/
│               ├── organic/
│               │   ├── extractors/
│               │   ├── processors/
│               │   ├── generators/
│               │   └── utils/
│               ├── inorganic/
│               │   ├── extractors/
│               │   ├── processors/
│               │   ├── generators/
│               │   └── utils/
│               └── ... (other subdomains)
├── data/
│   └── outputs/
│       └── chemistry/
│           ├── entity_dictionary.json
│           ├── graph_node_templates.json
│           ├── embedding_dictionary.json
│           └── subdomains/
│               ├── organic/
│               │   ├── entity_dictionary.json
│               │   ├── graph_node_templates.json
│               │   └── embedding_dictionary.json
│               ├── inorganic/
│               │   ├── entity_dictionary.json
│               │   ├── graph_node_templates.json
│               │   └── embedding_dictionary.json
│               └── ... (other subdomains)
├── interface/
│   └── streamlit_app.py
├── documentation/
│   └── HLD.md
└── requirements.txt                  # Updated with new dependencies
```

---

## 14. Usage Examples

### 14.1 Environment Setup

```python
from src.config.env_loader import EnvironmentLoader
from src.config.config_loader import ConfigLoader
from src.domains.chemistry.sources.source_factory import SourceFactory

# Setup environment (automatically detects or can be specified)
env_loader = EnvironmentLoader(domain="chemistry", environment="development")

# Create enhanced config loader
config_loader = ConfigLoader(domain="chemistry", environment="development")

# Create source factory with environment support
factory = SourceFactory(env_loader)

# Check available sources
available_sources = factory.get_available_sources()
print(f"Available sources: {available_sources}")

# Create PubChem source (if enabled in environment)
if env_loader.get('ENABLE_PUBCHEM', True):
    pubchem_source = factory.create_source_by_name('pubchem')
    print(f"PubChem settings: {env_loader.get_pubchem_settings()}")
```

### 14.2 Environment-Specific Configuration

```python
# Environment-specific settings are automatically loaded
if env_loader.is_development():
    print("Running in development mode with debug logging")
elif env_loader.is_production():
    print("Running in production mode with strict validation")
elif env_loader.is_testing():
    print("Running in testing mode with mock data sources")

# Access environment-specific configuration
storage_paths = config_loader.get_storage_paths()
rate_limiting = config_loader.get_rate_limiting_settings()
database_settings = config_loader.get_database_settings()
```

---

## 15. Downstream Use

- The final `.json` output (with provenance, review, UUIDs, and all occurrences) is the single source of truth for downstream extraction from publications and journals, ensuring traceability, quality, and human validation.
- UUIDs are assigned to each unique, reconciled entity instance only after extraction and reconciliation, not in the schema config.
- **Environment configuration** ensures consistent deployment across different environments with appropriate settings for each.

---

## 16. Next Steps

1. Approve and download this updated design document.
2. The environment configuration system is now fully implemented and ready for use.
3. Test the system with different environments (development, production, testing, staging).
4. Proceed to generate the folder structure, YAML configs, code modules, and example outputs for the chemistry domain.

---

**End of Document** 