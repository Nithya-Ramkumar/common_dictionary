# Detailed Design: Configuration Reconciliation Module

## Overview

The Configuration Reconciliation Module is a critical component of the Common Dictionary system that ensures all configuration files are consistent, complete, and properly cross-referenced before the system can proceed with data extraction and processing.

## Architecture

### Core Components

1. **ConfigReconciliation Class** - Main validation engine
2. **EnvironmentLoader Integration** - Flexible configuration path management
3. **Streamlit Interface** - User-friendly testing and validation interface
4. **Comprehensive Validation Logic** - Multi-level validation checks

### Design Principles

- **Environment-Aware**: Automatically adapts to different deployment environments
- **Flexible Input**: Supports both environment defaults and custom file uploads
- **Comprehensive Validation**: Multiple validation levels (syntax, structure, cross-references)
- **User-Friendly**: Clear error reporting and interactive testing interface
- **Extensible**: Easy to add new validation rules and configuration types

## Implementation Details

### 1. ConfigReconciliation Class

#### Key Methods

```python
class ConfigReconciliation:
    def __init__(self, env_loader=None, config_paths=None)
    def validate_configs(self) -> Dict[str, Any]
    def validate_entity_config(self, entity_config: Dict[str, Any]) -> bool
    def validate_source_mapping(self, source_mapping: Dict[str, Any]) -> bool
    def validate_conflict_resolution(self, conflict_resolution: Dict[str, Any]) -> bool
    def validate_cross_references(self, entity_attrs, source_priority_attrs, conflict_resolution_attrs) -> bool
    def get_validation_summary(self) -> str
```

#### Configuration Path Management

The system supports two modes of operation:

1. **Environment Default Mode**:
   - Uses EnvironmentLoader to get config paths from environment variables
   - Automatically adapts to different environments (dev, prod, test, staging)
   - Default behavior for normal operation

2. **Custom Path Mode**:
   - Allows overriding config paths for testing
   - Supports file uploads through Streamlit interface
   - Useful for testing with incomplete or modified configurations

#### Validation Levels

1. **File-Level Validation**:
   - File existence checks
   - YAML syntax validation
   - File content validation (non-empty, valid structure)

2. **Structure Validation**:
   - Entity configuration structure validation
   - Source mapping structure validation
   - Conflict resolution structure validation
   - Required field validation

3. **Cross-Reference Validation**:
   - Entity attributes vs source priority mappings
   - Entity attributes vs conflict resolution strategies
   - Orphaned mappings detection
   - Missing mappings detection

### 2. Environment Integration

#### EnvironmentLoader Usage

```python
# Initialize with environment
env_loader = EnvironmentLoader(environment="development")
reconciler = ConfigReconciliation(env_loader)

# Get config paths from environment
config_paths = {
    'entity_config': env_loader.get('ENTITY_CONFIG'),
    'source_mapping': env_loader.get('SOURCE_MAPPING'),
    'conflict_resolution': env_loader.get('CONFLICT_RESOLUTION'),
    'validation_config': env_loader.get('VALIDATION_CONFIG'),
    'extraction_config': env_loader.get('EXTRACTION_CONFIG'),
}
```

#### Environment-Specific Behavior

- **Development**: Detailed error messages, debug logging
- **Production**: Concise error reporting, performance optimized
- **Testing**: Strict validation, mock data support
- **Staging**: Intermediate settings for pre-production testing

### 3. Streamlit Interface

#### Features

1. **Environment Selection**:
   - Dropdown for environment selection (dev, prod, test, staging)
   - Real-time environment info display
   - Config path status indicators

2. **Configuration Source Options**:
   - Environment Default: Use existing config files
   - Upload Custom Files: Test with modified configurations

3. **File Upload Support**:
   - Individual YAML file uploads
   - Temporary file handling
   - Mixed mode (upload some, use defaults for others)

4. **Results Display**:
   - Tabbed interface (Summary, Errors, Warnings, Information)
   - Color-coded status indicators
   - Detailed error messages with context
   - Configuration file viewer

5. **Quick Actions**:
   - Test with incomplete configs
   - Export results to JSON
   - View current configuration files
   - Help and documentation

#### User Interface Layout

```
┌─────────────────────────────────────────────────────────────┐
│                    Configuration Validation                 │
├─────────────────────────────────────────────────────────────┤
│ Sidebar:                    │ Main Content:                │
│ - Environment Selection     │ - Validation Controls        │
│ - Config Source Options     │ - Results Display            │
│ - Environment Info          │ - File Upload (if custom)    │
│ - Config Path Status        │                              │
│ - Quick Actions             │                              │
│ - Help                      │                              │
└─────────────────────────────────────────────────────────────┘
```

## Validation Logic

### 1. Entity Configuration Validation

```python
def validate_entity_config(self, entity_config: Dict[str, Any]) -> bool:
    # Check for empty config
    # Validate entity structure
    # Check entity names
    # Validate attributes
    # Check required fields
    # Validate validation rules
```

**Validation Rules**:
- Config must not be empty
- Entities must have names
- Attributes must have names
- Required attributes must have validation rules
- Entity structure must be valid

### 2. Source Mapping Validation

```python
def validate_source_mapping(self, source_mapping: Dict[str, Any]) -> bool:
    # Check for empty mapping
    # Validate entity mappings
    # Check source configurations
    # Validate source lists
```

**Validation Rules**:
- Mapping must not be empty
- Entity types must be valid
- Source configurations must be complete
- Source lists must not be empty

### 3. Conflict Resolution Validation

```python
def validate_conflict_resolution(self, conflict_resolution: Dict[str, Any]) -> bool:
    # Check for empty resolution
    # Validate entity strategies
    # Check strategy configurations
```

**Validation Rules**:
- Resolution must not be empty
- Entity types must be valid
- Strategies must be properly configured

### 4. Cross-Reference Validation

```python
def validate_cross_references(self, entity_attrs, source_priority_attrs, conflict_resolution_attrs) -> bool:
    # Check missing mappings in source priority
    # Check missing mappings in conflict resolution
    # Check orphaned mappings in source priority
    # Check orphaned mappings in conflict resolution
```

**Validation Rules**:
- All entity attributes must have source priority mappings
- All entity attributes must have conflict resolution strategies
- No orphaned mappings in source priority
- No orphaned mappings in conflict resolution

## Error Handling and Reporting

### Error Categories

1. **Errors** (Critical):
   - File not found
   - YAML syntax errors
   - Missing required mappings
   - Orphaned mappings
   - Invalid structure

2. **Warnings** (Non-critical):
   - Empty source lists
   - Missing validation rules for required fields
   - Empty configurations

3. **Information** (Status):
   - Number of entities found
   - Number of attributes per entity
   - Validation statistics

### Error Reporting Format

```python
result = {
    'success': bool,
    'errors': List[str],
    'warnings': List[str],
    'info': List[str],
    'config_paths': Dict[str, str],
    'environment': str
}
```

### Error Message Examples

```
Errors:
1. [Source Priority] Entity 'Compound': Missing mappings for ['molecular_weight', 'boiling_point']
2. [Conflict Resolution] Entity 'Reaction': Missing mappings for ['catalyst']
3. File not found: /path/to/entity_config.yaml
4. YAML syntax error in source_mapping.yaml: line 5, column 10

Warnings:
1. Entity 'Compound' has no attributes
2. Required attribute 'name' in entity 'Compound' has no validation rules
3. Empty source list for 'formula' in entity 'Compound'

Information:
1. Found 3 entity types
2. - Compound: 5 attributes
3. - Reaction: 3 attributes
4. - Catalyst: 4 attributes
```

## Usage Examples

### 1. Command Line Usage

```bash
# Basic validation with default environment
python -m src.config_reconciliation

# Validation with specific environment
python -m src.config_reconciliation --env production

# Validation with custom config directory
python -m src.config_reconciliation --config-dir /custom/path

# Export results to JSON
python -m src.config_reconciliation --output results.json
```

### 2. Programmatic Usage

```python
from src.config.env_loader import EnvironmentLoader
from src.config.config_reconciliation import ConfigReconciliation

# Initialize with environment
env_loader = EnvironmentLoader(environment="development")
reconciler = ConfigReconciliation(env_loader)

# Run validation
result = reconciler.validate_configs()

# Check results
if result['success']:
    print("✅ All configurations are valid!")
else:
    print(f"❌ Found {len(result['errors'])} errors")
    for error in result['errors']:
        print(f"  - {error}")

# Get formatted summary
summary = reconciler.get_validation_summary()
print(summary)
```

### 3. Streamlit Interface Usage

```bash
# Run Streamlit interface
cd common_dictionary
streamlit run interface/streamlit_app.py
```

**Interface Features**:
- Select environment from dropdown
- Choose between environment defaults or custom file uploads
- Upload individual YAML files for testing
- View real-time validation results
- Export results to JSON
- View current configuration files

## Testing Strategy

### 1. Test Data Sets

**Valid Configurations** (Current configs):
- `common_dictionary/config/domains/chemistry/test_config/`
  - `entity_config.yaml` ✅
  - `source_priority.yaml` ✅
  - `conflict_resolution.yaml` ✅
  - `validation_config.yaml` ✅
  - `extraction_config.yaml` ✅

**Incomplete Configurations** (Test cases):
- `common_dictionary/tests/data/incomplete_configs/`
  - `entity_config_incomplete.yaml` (missing attributes)
  - `source_priority_incomplete.yaml` (orphaned mappings)
  - `conflict_resolution_incomplete.yaml` (missing strategies)

### 2. Test Scenarios

1. **Valid Configuration Test**:
   - Use current config files
   - Should pass all validations
   - Verify no errors or warnings

2. **Incomplete Configuration Test**:
   - Use intentionally broken configs
   - Should catch specific errors
   - Verify error messages are clear

3. **Environment Testing**:
   - Test with different environments
   - Verify environment-specific behavior
   - Check config path resolution

4. **File Upload Testing**:
   - Test custom file uploads
   - Verify mixed mode (upload + defaults)
   - Check temporary file handling

### 3. Validation Test Cases

**Entity Config Tests**:
- Empty config file
- Missing entity names
- Missing attribute names
- Required attributes without validation
- Invalid YAML syntax

**Source Mapping Tests**:
- Empty mapping file
- Missing entity types
- Empty source lists
- Invalid source configurations

**Conflict Resolution Tests**:
- Empty resolution file
- Missing entity types
- Missing strategies
- Invalid strategy configurations

**Cross-Reference Tests**:
- Missing mappings in source priority
- Missing mappings in conflict resolution
- Orphaned mappings in source priority
- Orphaned mappings in conflict resolution

## Integration Points

### 1. Environment System Integration

- Uses EnvironmentLoader for config path resolution
- Supports all environment types (dev, prod, test, staging)
- Integrates with environment-specific settings

### 2. Configuration System Integration

- Validates all configuration file types
- Ensures consistency across all config files
- Provides feedback for configuration improvements

### 3. Streamlit Interface Integration

- Provides user-friendly testing interface
- Supports interactive configuration testing
- Enables real-time validation feedback

### 4. Pipeline Integration

- Pre-validation step for data extraction pipeline
- Ensures system can proceed safely
- Provides detailed error reporting for debugging

## Future Enhancements

### 1. Additional Validation Rules

- Schema validation against JSON schemas
- Unit validation using Pint library
- Semantic validation of configuration values
- Best practices validation

### 2. Enhanced Error Reporting

- Error categorization by severity
- Suggested fixes for common errors
- Links to documentation
- Interactive error resolution

### 3. Configuration Generation

- Auto-generate missing configurations
- Suggest optimal configurations
- Configuration templates
- Migration tools for config updates

### 4. Performance Optimizations

- Parallel validation for large configs
- Caching of validation results
- Incremental validation
- Background validation

## Conclusion

The Configuration Reconciliation Module provides a robust foundation for ensuring configuration consistency and completeness in the Common Dictionary system. Its environment-aware design, comprehensive validation logic, and user-friendly interface make it an essential component for maintaining system reliability and facilitating development and testing workflows.

The module's extensible architecture allows for future enhancements while maintaining backward compatibility and providing clear, actionable feedback for configuration issues. 