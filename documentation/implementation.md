# Implementation & Testing Plan: Common Dictionary Module

---

## Progress Table

| Step/Module           | Status   |
|-----------------------|----------|
| config_reconciliation | ⬜ To Do  |
| extractor             | ⬜ To Do  |
| reconciler            | ⬜ To Do  |
| uuid_assigner         | ⬜ To Do  |
| validator             | ⬜ To Do  |
| enricher              | ⬜ To Do  |
| synonym_clusterer     | ⬜ To Do  |
| review_logger         | ⬜ To Do  |
| exporter              | ⬜ To Do  |
| unit_converter        | ⬜ To Do  |
| Streamlit UI          | ⬜ To Do  |
| End-to-End            | ⬜ To Do  |

---

## Step-by-Step Implementation & Testing Plan

### 1. Initial Setup

**A. Directory Structure**
- Ensure your folders are set up as per the design (`config/`, `data/`, `src/`, `interface/`, `documentation/`).

**B. Environment**
- Set up a Python virtual environment.
- Install core dependencies:  
  `pandas`, `PyYAML`, `jsonschema`, `cerberus`, `pint`, `scikit-learn`, `networkx`, `sentence-transformers`, `transformers`, `spacy`, `streamlit`, `pyvis`, `umap-learn`, `plotly`, etc.

**Test/Verify:**  
- All dependencies install without error.
- Directory structure matches the design.

---

### 2. Module-by-Module Implementation

#### 2.1. config_reconciliation.py
- **Purpose:** Validate that all entity attributes in `entity_config.yaml` are mapped in `source_priority.yaml` and `conflict_resolution.yaml`.
- **Test/Verify:**
  - Run with intentionally incomplete configs and confirm errors are caught.
  - Run with correct configs and confirm a clean report.

#### 2.2. extractor.py
- **Purpose:** Extract raw entities and relationships from sources, with provenance.
- **Test/Verify:**
  - Output files are generated and match expected schema.
  - Provenance is included for each value.
  - Test with both dummy and (eventually) real data.

#### 2.3. reconciler.py
- **Purpose:** Merge, deduplicate, and reconcile extracted data using source priority and conflict rules.
- **Test/Verify:**
  - Conflicting values are resolved as per config.
  - All alternatives and provenance are preserved.
  - Edge cases: multiple sources, missing data, etc.

#### 2.4. uuid_assigner.py
- **Purpose:** Assign UUIDs to unique entities and relationships.
- **Test/Verify:**
  - All entities/relationships have unique, valid UUIDs.
  - No duplicates.

#### 2.5. validator.py
- **Purpose:** Validate data against schema, units, value ranges, synonyms, etc.
- **Test/Verify:**
  - Invalid data is flagged with clear errors.
  - Valid data passes all checks.
  - Test with edge cases (bad units, out-of-range values, etc.).

#### 2.6. enricher.py
- **Purpose:** Enrich entities/relationships with additional data from APIs or local tools.
- **Test/Verify:**
  - Enrichment fields are added as specified.
  - Provenance for enrichment is tracked.

#### 2.7. synonym_clusterer.py
- **Purpose:** Cluster and validate synonyms for each entity.
- **Test/Verify:**
  - Synonyms are clustered as expected.
  - Canonical terms are selected per config.
  - Test with known synonym sets.

#### 2.8. review_logger.py
- **Purpose:** Log all human-in-the-loop review actions from the UI.
- **Test/Verify:**
  - All review actions are logged with timestamp, user, and action.
  - Review history is correctly attached to entities/relationships.

#### 2.9. exporter.py
- **Purpose:** Export final, reviewed, validated data for downstream use.
- **Test/Verify:**
  - Output files match the documented schema.
  - All provenance, review history, and UUIDs are present.

#### 2.10. unit_converter.py (Optional/Advanced)
- **Purpose:** Validate and convert units using Pint.
- **Test/Verify:**
  - Units are converted to canonical forms.
  - Invalid units are flagged.

---

### 3. Streamlit Interface
- **Purpose:** Modular UI for each module, with output display, config editing, rerun, visualization, and review logging.
- **Test/Verify:**
  - UI loads and displays data.
  - Editing and rerunning works as expected.
  - Visualizations are correct and interactive.
  - Review actions are logged and visible.

---

### 4. End-to-End Testing
- **Purpose:** Validate the full pipeline with a small, known dataset.
- **Test/Verify:**
  - Each output matches expectations at every stage.
  - Provenance, review history, and all metadata are preserved.
  - Error handling and edge cases are covered.

---

## Summary Table: Implementation & Testing

| Step/Module           | What to Implement                | What to Test/Verify                                 |
|-----------------------|----------------------------------|-----------------------------------------------------|
| config_reconciliation | Config validation logic          | Catches missing/orphaned mappings                   |
| extractor             | Extraction logic, provenance     | Output schema, provenance, dummy/real data          |
| reconciler            | Merge, dedup, conflict resolve   | Correct resolution, alternatives, provenance        |
| uuid_assigner         | UUID assignment                  | All entities/relations have unique UUIDs            |
| validator             | Schema/unit/synonym validation   | Flags invalid, passes valid, edge cases             |
| enricher              | Enrichment logic                 | Fields added, provenance tracked                    |
| synonym_clusterer     | Embedding, clustering            | Correct clusters, canonical terms                   |
| review_logger         | Logging, UI integration          | All actions logged, attached to data                |
| exporter              | Output generation                | Schema match, all metadata present                  |
| unit_converter        | Unit conversion                  | Canonical units, invalid flagged                    |
| Streamlit UI          | Modular UI, rerun, edit, viz     | UI works, rerun/edit/viz/review all functional      |
| End-to-End            | Full pipeline                    | All outputs, metadata, error handling, edge cases   | 