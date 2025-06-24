# Common Dictionary Module: Detailed Design Document

---

## 1. Overview

The Common Dictionary Module is a modular, extensible system for extracting, reconciling, validating, enriching, and reviewing entities and relationships in scientific domains (starting with chemistry). It emphasizes provenance, traceability, human-in-the-loop review, and standards-compliant outputs for downstream applications such as knowledge graph construction and embedding-based search.

---

## 2. Architecture & Module Pipeline

### 2.1. Pipeline Modules

| Step | Module                  | Purpose                                                                                  | Input Files/Source                | Output Files/Artifacts                |
|------|-------------------------|-----------------------------------------------------------------------------------------|-----------------------------------|---------------------------------------|
| 1    | config_reconciliation.py| Validate config consistency and completeness                                            | entity_config.yaml, source_priority.yaml, conflict_resolution.yaml | Validation report, logs               |
| 2    | extractor.py            | Extract raw entities and relationships from sources, with provenance                    | extraction_config.yaml, source APIs, relationship.yaml (template) | raw_extracted_data.json, raw_relationships.json |
| 3    | reconciler.py           | Merge, deduplicate, and reconcile extracted data using source priority and conflict rules| raw_extracted_data.json, raw_relationships.json, source_priority.yaml, conflict_resolution.yaml | reconciled_entities.json, reconciled_relationships.json |
| 4    | uuid_assigner.py        | Assign UUIDs to unique entities and relationships                                       | reconciled_entities.json, reconciled_relationships.json | entities_with_uuids.json, relationships_with_uuids.json |
| 5    | validator.py            | Validate data against schema, units, value ranges, synonyms, etc.                       | validation_config.yaml, metric_units.yaml, entities_with_uuids.json, relationships_with_uuids.json | Validation report, cleaned data        |
| 6    | enricher.py             | Enrich entities/relationships with additional data from APIs or local tools             | entities_with_uuids.json, relationships_with_uuids.json, enrichment rules | enriched_entities.json, enriched_relationships.json |
| 7    | synonym_clusterer.py    | Cluster and validate synonyms for each entity                                           | entity data, validation_config.yaml | synonym_clusters.json                  |
| 8    | review_logger.py        | Log all human-in-the-loop review actions from the UI                                    | Review actions (from UI)           | review_history.json, updated data      |
| 9    | exporter.py             | Export final, reviewed, validated data for downstream use (entities, relationships, graphs, embeddings) | Final entity/relationship data     | entity_dictionary.json, relationship_dictionary.json, graph.json, embeddings.json |
| 10   | unit_converter.py*      | (Optional) Runtime unit validation and conversion using Pint                            | Data with units, metric_units.yaml | Data with canonicalized units          |

*Optional/Advanced

---

### 2.2. Data & Config Structure

- **entity_config.yaml**: Defines entity types, attributes, validation, enrichment rules.
- **relationship.yaml**: Template for relationship schema (actual relationships are extracted).
- **extraction_config.yaml**: Source connection details, extraction parameters.
- **source_priority.yaml**: Source priority for each attribute.
- **conflict_resolution.yaml**: Conflict resolution strategies.
- **validation_config.yaml**: Validation rules (schema, units, synonyms, etc.).
- **metric_units.yaml**: SI units, allowed units, conversion rules.

---

### 2.3. Output Structure

- **entities_with_uuids.json**:  
  Each entity references its class/type (e.g., `"type": "Compound"` for water), includes UUID, attributes, provenance, review history, and relationships (with target UUID/type).

- **relationships_with_uuids.json**:  
  Each relationship includes UUID, type, from/to entity UUIDs, attributes, provenance, review history.

- **graph.json**:  
  Nodes (entities) and edges (relationships) with all metadata, for downstream graph construction.

- **embeddings.json**:  
  Embedding vectors for each node/edge, with model info and provenance.

---

## 3. Streamlit Interface Design

### 3.1. Purpose

- For each module, provide a UI to:
  - View output (entities, relationships, graphs, embeddings)
  - Edit input/config (YAML/JSON, parameters)
  - Rerun the module
  - Visualize graphs and embeddings
  - Log and display review actions

### 3.2. Structure

- `interface/streamlit_app.py`: Main entry point.
- `interface/modules/`: One file per module (e.g., `extraction.py`, `validation.py`, `graph.py`).
- **Sidebar**: Module/domain selection.
- **Main Area**: Output display, config editing, rerun button, visualization, review history.

### 3.3. Features

- **Config Editor**: YAML/JSON editor with validation.
- **Output Viewer**: JSON/tables for entities/relationships, interactive graph and embedding visualizations.
- **Rerun Button**: Triggers backend module execution.
- **Review Logging**: All actions logged and displayed.
- **Download/Upload**: For configs, outputs, embeddings.
- **Visualization**:  
  - Graphs: `networkx`, `pyvis`, or `streamlit-agraph`
  - Embeddings: `plotly`, `matplotlib`, `umap-learn`, `scikit-learn` (PCA/TSNE)

---

## 4. LLMs and Advanced AI Tools

### 4.1. Use Cases for LLMs

- **Entity/Relationship Extraction**: Use locally hosted LLMs (e.g., HuggingFace Transformers, spaCy) to extract entities and relationships from unstructured text, supplementing rule-based or API-based extraction.
- **Synonym Detection/Clustering**: Use local LLM embeddings (e.g., Sentence Transformers) for semantic similarity in synonym clustering.
- **Enrichment**: Use local LLMs to generate descriptions, summaries, or infer missing attributes for entities/relationships.
- **Validation/Reasoning**: Use local LLMs for cross-field logic checks, anomaly detection, or suggesting corrections during human review.

### 4.2. Integration Points

- **extractor.py**: Optionally call local LLM models for text extraction.
- **synonym_clusterer.py**: Use local LLM-based embeddings for clustering.
- **enricher.py**: Use local LLMs for enrichment tasks.
- **validator.py**: Use local LLMs for advanced validation or reasoning.

### 4.3. LLM/AI Libraries (Local Only)

- **HuggingFace Transformers**: For running local LLMs and embedding models ([transformers](https://huggingface.co/docs/transformers/))
- **Sentence Transformers**: For semantic embeddings ([sentence-transformers](https://www.sbert.net/))
- **spaCy**: For NER and linguistic features ([spaCy](https://spacy.io/))
- **llama.cpp, vLLM, or similar**: For efficient local LLM inference ([llama.cpp](https://github.com/ggerganov/llama.cpp), [vLLM](https://github.com/vllm-project/vllm))

> **Note:** All LLMs and embedding models must be hosted and run locally. No cloud APIs (e.g., OpenAI API) are to be used. Choose models that fit your hardware and privacy requirements. Use quantized models or optimized inference engines (e.g., llama.cpp, vLLM) for efficiency if needed.

---

## 5. Tools and Libraries

### 5.1. Core Python Libraries

| Purpose                        | Library/Tool         | Notes/Links                                      |
|--------------------------------|----------------------|--------------------------------------------------|
| Config parsing (YAML/JSON)     | PyYAML, json         | https://pyyaml.org/                              |
| Data manipulation              | pandas               | https://pandas.pydata.org/                       |
| UUID generation                | uuid                 | Standard library                                 |
| Validation (schema, regex)     | Cerberus, jsonschema | https://docs.python-cerberus.org/en/stable/      |
| Unit handling/conversion       | pint                 | https://pint.readthedocs.io/                     |
| API requests                   | requests, httpx      | https://docs.python-requests.org/                |
| Provenance tracking            | Custom, or `prov`    | https://pypi.org/project/prov/                   |
| Clustering (synonyms)          | scikit-learn, numpy  | https://scikit-learn.org/                        |
| Graph construction             | networkx             | https://networkx.org/                            |
| Embedding generation           | sentence-transformers, transformers | https://www.sbert.net/           |
| Visualization (graphs)         | pyvis, streamlit-agraph | https://pyvis.readthedocs.io/                |
| Visualization (embeddings)     | plotly, matplotlib, umap-learn | https://plotly.com/python/umap/           |
| Streamlit UI                   | streamlit            | https://streamlit.io/                            |
| Logging                        | logging              | Standard library                                 |
| Human review DB (optional)     | sqlite3, TinyDB      | For persistent review logs                       |
| LLM/AI Libraries (Local Only)  | transformers, sentence-transformers, spaCy, llama.cpp, vLLM | See above section |

---

### 5.2. Table: Tools and Libraries

| Module/Functionality           | Library/Tool           | Purpose/Usage                                    |
|-------------------------------|------------------------|--------------------------------------------------|
| Config parsing                 | PyYAML, json           | Read/write YAML/JSON configs                     |
| Data manipulation              | pandas                 | Tabular data processing                          |
| UUID assignment                | uuid                   | Generate unique IDs                              |
| Schema/validation              | Cerberus, jsonschema   | Validate data against schema/rules                |
| Unit handling                  | pint                   | SI units, conversions                            |
| API access                     | requests, httpx        | Fetch data from external sources                 |
| Provenance tracking            | prov, custom           | Track and store provenance metadata              |
| Clustering (synonyms)          | scikit-learn, numpy    | Synonym clustering (e.g., DBSCAN)                |
| Graph construction             | networkx               | Build and manipulate graphs                      |
| Embedding generation           | sentence-transformers, transformers | Generate embeddings for entities/relations |
| Graph visualization            | pyvis, streamlit-agraph| Interactive graph display in UI                  |
| Embedding visualization        | plotly, matplotlib, umap-learn | Visualize high-dimensional embeddings    |
| Streamlit UI                   | streamlit              | Build modular, interactive UI                    |
| Logging                        | logging                | Log actions, errors, reviews                     |
| Review DB (optional)           | sqlite3, TinyDB        | Store review history                             |
| LLM/AI Libraries (Local Only)  | transformers, sentence-transformers, spaCy, llama.cpp, vLLM | Local LLM-based extraction, enrichment, validation |

---

## 6. Example Output Structure

### 6.1. Entity Output Example
```json
{
  "uuid": "1234-5678-...",
  "type": "Compound",
  "name": "Water",
  "attributes": { ... },
  "provenance": { ... },
  "review_history": [ ... ],
  "relationships": [
    {
      "type": "participates_in",
      "target_uuid": "5678-1234-...",
      "target_type": "Reaction"
    }
  ]
}
```

### 6.2. Relationship Output Example
```json
{
  "uuid": "rel-9876-...",
  "type": "catalyzes",
  "from_uuid": "cat-1234-...",
  "from_type": "Catalyst",
  "to_uuid": "rxn-5678-...",
  "to_type": "Reaction",
  "attributes": { ... },
  "provenance": { ... },
  "review_history": [ ... ]
}
```

### 6.3. Graph Output Example
```json
{
  "nodes": [ ... ],  // Entities with UUIDs, types, attributes
  "edges": [ ... ]   // Relationships with UUIDs, types, from/to, attributes
}
```

### 6.4. Embeddings Output Example
```json
{
  "entities": {
    "1234-5678-...": {
      "embedding": [0.1, 0.2, ...],
      "model": "chem-bert-v1",
      "provenance": { ... }
    }
  },
  "relationships": {
    "rel-9876-...": {
      "embedding": [0.3, 0.4, ...],
      "model": "rel-bert-v1",
      "provenance": { ... }
    }
  }
}
```

---

## 7. Extensibility & Best Practices

- **Domain Extensibility:**  
  Add new domains/subdomains by copying config/output/module structure.
- **Traceability:**  
  All edits, reviews, and reruns are logged and versioned.
- **Validation:**  
  All config edits are schema-validated before rerun.
- **Visualization:**  
  Use interactive visualizations for graphs and embeddings to aid review and exploration.
- **Provenance:**  
  All data points, relationships, and embeddings include provenance and review history.
- **LLM/AI Integration:**  
  LLMs and advanced AI tools can be plugged in at extraction, enrichment, synonym clustering, and validation stages for improved automation and intelligence.

---

## 8. Next Steps

- Confirm or adjust the module list, output structure, and tool/library choices.
- Prioritize which scripts/modules to implement first.
- Optionally, generate code templates for foundational modules (e.g., `config_reconciliation.py`, `extractor.py`, Streamlit UI skeleton).

---

**This design ensures modularity, traceability, extensibility, and a robust human-in-the-loop workflow, with all outputs ready for downstream graph and embedding applications, and is future-proofed for LLM/AI integration.** 