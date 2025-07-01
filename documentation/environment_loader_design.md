# EnvironmentLoader Design and Usage

## Purpose and Design Principles
- **Centralized loading**: Loads all environment variables for the chemistry domain (or any specified domain/environment) from a single `.env` file.
- **Single source of truth**: Ensures all modules access the same environment configuration.
- **Flexibility**: Loads all variables, not just a hardcoded subset, so new variables can be added without changing the loader.
- **Portability**: Works from any directory, as long as absolute paths are used in the `.env` file for file-based configs.

## How It Works
- On initialization, loads the appropriate `.env` file (e.g., `env.development`) from the `env_templates/<domain>/` directory.
- Loads all variables from the environment into an internal dictionary (`self.env_vars`).
- Provides methods to get a single variable (`get(key)`), get all variables (`get_all()`), and print all variables (`print_config_paths()`).

## Usage in Code
```python
from config.env_loader import EnvironmentLoader

# Initialize for the chemistry domain (default: development)
env_loader = EnvironmentLoader(domain="chemistry")

# Get a variable
entity_config_path = env_loader.get('ENTITY_CONFIG')

# Get all variables
env_vars = env_loader.get_all()

# Print all variables (for debugging)
env_loader.print_config_paths()
```

## Best Practices for .env Files
- Use **absolute paths** for all file-based config variables (e.g., `ENTITY_CONFIG`, `ONTOLOGY_CONFIG`, etc.) to ensure portability.
- Keep all environment variables for a given environment in a single `.env` file (per environment).
- Document any required variables in your README or in comments at the top of the `.env` file.

## Example .env Snippet
```env
ENTITY_CONFIG=/abs/path/to/entity_config.yaml
ONTOLOGY_CONFIG=/abs/path/to/ontology.yaml
METRIC_UNITS=/abs/path/to/metric_units.yaml
DB_HOST=localhost
DB_PORT=5432
# ... other variables ...
```

## Extensibility Notes
- You can add new variables to your `.env` file at any time; they will be available via the loader without code changes.
- If you need to support multiple domains or environments, instantiate the loader with the appropriate arguments.
- For large projects, you can add helper methods to group or validate related variables. 