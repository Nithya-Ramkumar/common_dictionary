# AI Coding Instructions for This Directory

When making any changes in this directory, **always follow these instructions**:

1. **Review Config Files**  
    - Always check config files, especially entity and source mapping in `confign/domains/chemistry/testconfig2` and other latest directories.

2. **Config File Changes**  
    - When suggesting changes to config files, consider the impact on `config_reconciliation.py` and any related feature updates.

3. **Source Files Management**  
    - Source files for each extraction source are in `src/domain/chemistry/sources/`.
    - When creating or editing these files, ensure there is a mapping section from the entity field to the corresponding source files.
    - All source files must align with `base_source.py` and consistently implement the `search` and `extract_by_key` methods.

4. **Output Generator Updates**  
    - If you update `outputgenerator.py`, also update `env_templates/config_templates/output_entity_schema.template.yaml` with the generic schema for downstream interpretation (cross-domain).

5. **Entity File Updates**  
    - When updating entity files, also update the template files in `env_templates/config_templates/`.
    - These templates are generic and cross-domain; add/update them with clear comments to help users create specific files.

6. **Debugging Standards**  
    - Always review the `env.testing` file for debug flags.
    - All debugging messages must follow the specified format and include relevant debug information.

**Follow these steps to ensure consistency, maintainability, and cross-domain compatibility.**