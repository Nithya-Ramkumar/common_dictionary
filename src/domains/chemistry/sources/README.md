# Chemistry Sources Directory

This directory contains all source extractor classes for the chemistry domain. Each source is responsible for extracting entities and relationships from a specific data source (API, database, file, local library, etc.).

## Adding a New Source
1. Create a new Python file (e.g., `mydb_source.py`) and implement a class inheriting from `BaseSource`.
2. At the end of your source file, register your class with the source factory:
   ```python
   from .source_factory import SourceFactory
   SourceFactory.register_source_type("mydb", MyDBSource)
   ```
3. Add your source to the `extraction_config.yaml` under `extraction_sources` with `type: mydb` and the appropriate connection/schema details.

## Extraction Config Reference
- See `../../../../env_templates/chemistry/config_templates/extraction_config.template.yaml` for the full template and field descriptions.
- Each source in the config must have:
  - `name`: Unique name
  - `type`: Must match the string used in `register_source_type`
  - `enabled`: true/false
  - `description`: Short description
  - `connection`: Connection details (API URL, DB credentials, file path, etc.)
  - `schema`: Endpoints/tables/fields and their types

## Source Factory
- The factory uses a registry/plugin pattern to support unlimited source types.
- When the orchestrator loads the extraction config, it uses the factory to instantiate each enabled source.
- No changes to the factory are needed to add new sources—just register them!

## More Information
- See the detailed design documentation: `../../../../documentation/detailed_design_config_recon.md`
- See the extraction config template for examples and field descriptions. 