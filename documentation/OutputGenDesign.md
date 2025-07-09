# Output Generator Design Documentation

## Overview

The OutputGenerator is a core component of the chemistry data extraction pipeline that generates structured output files in multiple formats (JSON, CSV, TSV, HTML) from extracted entity data. It implements a **flat attribute record approach** that ensures per-attribute timestamp and provenance tracking while maintaining both human readability and machine processability.

## Key Design Principles

### 1. Flat Attribute Records
- **One row per attribute extraction**: Each extracted attribute becomes a separate record
- **Per-attribute metadata**: Each record includes its own timestamp, provenance, and confidence score
- **Traceability**: Complete audit trail for every data point
- **Flexibility**: Easy to filter, group, or analyze individual attributes

### 2. Dual Output Strategy
- **Machine-Processable Formats** (JSON/CSV/TSV): Flat attribute records for downstream processing
- **Human-Readable Format** (HTML): CID-grouped visual display with per-attribute metadata

## Output Format Specifications

### JSON Output
**Structure**: Flat list of attribute records with summary metadata

```json
{
  "summary": {
    "polymers_requested": 20,
    "polymers_extracted": 15,
    "timestamp": "2024-07-07 18:10:00",
    "total_entities": 75,
    "entity_type_counts": {"Polymer": 15},
    "source_counts": {"pubchem": 75}
  },
  "entities": [
    {
      "cid": "12345",
      "attribute": "name",
      "value": "Polymer X",
      "timestamp": "2024-07-07 18:10:00",
      "provenance": {
        "source": "pubchem",
        "endpoint": "property",
        "query": {"cid": "12345"}
      },
      "confidence": 1.0,
      "source": "pubchem"
    },
    {
      "cid": "12345",
      "attribute": "molecular_weight",
      "value": 1234.5,
      "timestamp": "2024-07-07 18:10:01",
      "provenance": {
        "source": "pubchem",
        "endpoint": "property",
        "query": {"cid": "12345"}
      },
      "confidence": 1.0,
      "source": "pubchem"
    }
  ]
}
```

**Pros**: 
- Best for programmatic processing
- Flexible structure
- Easy to load into pandas, SQL, or other tools

**Cons**: 
- Less compact for human browsing

### CSV/TSV Output
**Structure**: One row per attribute extraction with standardized columns

```csv
# polymers_requested: 20
# polymers_extracted: 15
# timestamp: 2024-07-07 18:10:00
cid,attribute,value,timestamp,provenance,confidence,source
12345,name,Polymer X,2024-07-07 18:10:00,pubchem | endpoint=property | query=(cid=12345),1.0,pubchem
12345,molecular_weight,1234.5,2024-07-07 18:10:01,pubchem | endpoint=property | query=(cid=12345),1.0,pubchem
```

**Columns**:
- `cid`: Compound ID from PubChem
- `attribute`: Attribute name (e.g., "name", "molecular_weight")
- `value`: Extracted value
- `timestamp`: Extraction timestamp for each attribute
- `provenance`: Formatted provenance string
- `confidence`: Extraction confidence score
- `source`: Data source name

**Pros**:
- Easy for Excel/pandas processing
- Clear structure
- Standard format for data analysis

**Cons**:
- More rows than entity-grouped approach
- Provenance information is string-formatted

### HTML Output
**Structure**: CID-grouped visual display with per-attribute metadata

```html
<!DOCTYPE html>
<html>
<head>
  <title>Polymer Extraction Report</title>
  <!-- CSS styling for tabs and tables -->
</head>
<body>
  <h1>Polymer Extraction Report</h1>
  <div><b>Run timestamp:</b> 2024-07-07 18:10:00</div>
  
  <h2>Summary</h2>
  <ul>
    <li><b>Polymers requested:</b> 20</li>
    <li><b>Polymers extracted:</b> 15</li>
    <li><b>Errors:</b> 0</li>
  </ul>
  
  <h2>Extracted Polymers (Grouped by CID)</h2>
  
  <!-- Navigation tabs -->
  <div class="nav-tabs">
    <div class="nav-tab active" onclick="showTab('cid-12345')">CID 12345</div>
    <div class="nav-tab" onclick="showTab('cid-67890')">CID 67890</div>
  </div>
  
  <!-- Tab content for each CID -->
  <div id="cid-12345" class="tab-content active">
    <div class="cid-section">
      <div class="cid-header"><h3>CID: 12345</h3></div>
      <table class="attribute-table">
        <tr><th>Attribute</th><th>Value</th><th>Source</th><th>Timestamp</th><th>Confidence</th></tr>
        <tr>
          <td><strong>name</strong></td>
          <td>Polymer X</td>
          <td class="provenance">pubchem: cid=12345</td>
          <td>2024-07-07 18:10:00</td>
          <td>1.0</td>
        </tr>
        <tr>
          <td><strong>molecular_weight</strong></td>
          <td>1234.5</td>
          <td class="provenance">pubchem: cid=12345</td>
          <td>2024-07-07 18:10:01</td>
          <td>1.0</td>
        </tr>
      </table>
    </div>
  </div>
</body>
</html>
```

**Features**:
- **Tabbed Navigation**: Each CID gets its own tab
- **Visual Grouping**: Clear separation between different compounds
- **Per-Attribute Metadata**: Each attribute shows source, timestamp, and confidence
- **Responsive Design**: Clean, modern styling
- **Interactive Tabs**: JavaScript functionality for navigation

**Pros**:
- Human-friendly and readable
- Easy to scan and review
- Visual organization by CID
- Complete metadata visibility

**Cons**:
- Large file size for big datasets
- Not suitable for programmatic processing

## Schema Configuration

### Output Entity Schema Structure

The output schema is defined in `output_entity_schema.yaml` and supports multiple schema types:

```yaml
output_schemas:
  default:
    - cid                    # Compound ID from PubChem
    - attribute              # Attribute name
    - value                  # Extracted value
    - timestamp              # Extraction timestamp for each attribute
    - provenance             # Source and extraction details
    - confidence             # Extraction confidence score
    - source                 # Data source name
  
  summary:
    - cid
    - attribute
    - value
    - timestamp
    - source
  
  detailed:
    - cid
    - attribute
    - value
    - timestamp
    - provenance
    - confidence
    - source
    - entity_type
    - validation_status
    - extraction_method
    - version
```

### Schema Types

1. **default**: Standard output with all key fields for most use cases
2. **summary**: Minimal output for quick overviews or dashboards
3. **detailed**: Comprehensive output for in-depth analysis or auditing

## Implementation Details

### OutputGenerator Class

```python
class OutputGenerator:
    """
    Generates output files (JSON, CSV, TSV, HTML) for extracted entities using a schema.
    
    Output formats:
    - JSON/CSV/TSV: Flat attribute records (one row per attribute extraction)
    - HTML: CID-grouped display with per-attribute metadata
    """
    
    def generate_outputs(self, entities, summary, output_dir, schema_name='default', formats=['json', 'csv', 'tsv']):
        """Generate outputs in specified formats using flat attribute records."""
        
    def _write_json(self, entities, summary, output_dir):
        """Write JSON output as flat attribute records with per-attribute metadata."""
        
    def _write_table(self, entities, summary, output_dir, fields, sep, ext):
        """Write CSV/TSV output as flat attribute records (one row per attribute extraction)."""
        
    def _write_html(self, entities, summary, output_dir, fields):
        """Write HTML output with CID-grouped display and per-attribute metadata."""
        
    def _format_provenance_for_csv(self, provenance):
        """Format provenance dict as a string for CSV output."""
```

### Key Methods

1. **`_write_json()`**: Outputs flat attribute records in JSON format
2. **`_write_table()`**: Outputs flat attribute records in CSV/TSV format with formatted provenance
3. **`_write_html()`**: Creates CID-grouped HTML display with per-attribute metadata
4. **`_format_provenance_for_csv()`**: Converts provenance dict to readable string format

## Data Flow

### Input Data Structure
```python
entities = [
    {
        'cid': '12345',
        'attribute': 'name',
        'value': 'Polymer X',
        'timestamp': '2024-07-07 18:10:00',
        'provenance': {'source': 'pubchem', 'endpoint': 'property', 'query': {'cid': '12345'}},
        'confidence': 1.0,
        'source': 'pubchem'
    },
    # ... more attribute records
]
```

### Processing Flow
1. **Input**: Flat attribute records from extractor
2. **JSON/CSV/TSV**: Output records as-is (flat structure)
3. **HTML**: Group by CID for visual display, but maintain per-attribute metadata

## Benefits of Flat Attribute Records

### 1. Traceability
- Every attribute has its own timestamp and provenance
- Complete audit trail for data extraction
- Easy to track when and where each piece of data came from

### 2. Flexibility
- Easy to filter by attribute, source, or confidence
- Simple to group or aggregate in downstream tools
- Supports complex queries and analysis

### 3. Machine Processability
- Natural fit for pandas DataFrames
- Easy to load into SQL databases
- Compatible with data analysis tools

### 4. Human Readability
- HTML output provides visual organization
- Clear separation of concerns
- Easy to scan and review

## Comparison with Entity-Grouped Approach

| Aspect | Flat Attribute Records | Entity-Grouped Records |
|--------|----------------------|------------------------|
| **Traceability** | ✅ Per-attribute metadata | ❌ Entity-level metadata only |
| **Processing** | ✅ Easy to filter/group | ❌ Requires data transformation |
| **Audit Trail** | ✅ Complete | ❌ Limited |
| **Human Reading** | ✅ HTML groups by CID | ✅ Natural entity view |
| **File Size** | ⚠️ More rows | ✅ Fewer rows |
| **Flexibility** | ✅ High | ❌ Low |

## Usage Examples

### Loading into Pandas
```python
import pandas as pd

# Load CSV output
df = pd.read_csv('entities.csv')
print(df.head())

# Filter by attribute
name_data = df[df['attribute'] == 'name']

# Group by CID
by_cid = df.groupby('cid')

# Filter by confidence
high_confidence = df[df['confidence'] > 0.8]
```

### Loading into SQL
```sql
-- Create table from CSV
CREATE TABLE extracted_attributes (
    cid VARCHAR(20),
    attribute VARCHAR(50),
    value TEXT,
    timestamp TIMESTAMP,
    provenance TEXT,
    confidence FLOAT,
    source VARCHAR(20)
);

-- Query examples
SELECT * FROM extracted_attributes WHERE attribute = 'name';
SELECT cid, COUNT(*) as attr_count FROM extracted_attributes GROUP BY cid;
SELECT * FROM extracted_attributes WHERE confidence > 0.8;
```

## Configuration Files

### Main Schema
- **Location**: `common_dictionary/config/domains/chemistry/output_entity_schema.yaml`
- **Purpose**: Defines output structure for main chemistry domain

### Test Config Schema
- **Location**: `common_dictionary/config/domains/chemistry/test_config1/output_entity_schema.yaml`
- **Purpose**: Defines output structure for test configurations

## Future Enhancements

### Potential Improvements
1. **Compression**: Add gzip compression for large outputs
2. **Streaming**: Support for streaming large datasets
3. **Custom Formats**: Add support for additional output formats (XML, Parquet, etc.)
4. **Validation**: Add output validation against schema
5. **Caching**: Implement caching for repeated queries

### Extensibility
- Schema-driven approach allows easy addition of new fields
- Modular design supports new output formats
- Configuration-based approach enables domain-specific customization

## Conclusion

The OutputGenerator implements a robust, flexible approach to data output that balances human readability with machine processability. The flat attribute record structure ensures complete traceability and auditability while maintaining ease of use for both technical and non-technical users.

The dual output strategy (flat records for processing, grouped display for reading) provides the best of both worlds, making the extracted data accessible and useful for a wide range of downstream applications. 