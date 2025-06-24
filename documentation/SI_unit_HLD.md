# SI Units, Multiples, and Derived Units — High-Level Design (HLD)

---

## 1. Overview

This document describes the design and approach for handling quantitative fields, SI units, multiples (prefixes), and derived units in the Common Dictionary Module and across the entire project. The goal is to ensure modularity, extensibility, and compliance with industry standards (SI, NIST, BIPM) for all quantitative data.

---

## 2. Centralized, Shared Units File

- A single `metric_units.yaml` file is generated from the Pint library (or another SI-compliant source) and stored outside the `common_dictionary` folder, e.g., in a `config/` directory at the project root.
- This file serves as the authoritative list of allowed SI and metric units, including multiples (prefixes) and derived units.
- All modules (not just the common dictionary) reference this file for unit validation and conversion, ensuring consistency project-wide.

---

## 3. Entity Design: Quantitative Fields and Units

- Entities with quantitative attributes (e.g., mass, length, concentration) have:
  - A value field (e.g., `mass_value`)
  - A unit field (e.g., `mass_unit`)
- Unit fields do not hardcode allowed values; they are validated against the central `metric_units.yaml`.

**Example (entity YAML):**
```yaml
attributes:
  - name: mass_value
    type: float
    required: false
  - name: mass_unit
    type: string
    required: false
```

---

## 4. Validation Config

- Validation rules for unit fields reference the central `metric_units.yaml` for allowed units.
- No SI units are hardcoded in the validation config; all validation is dynamic and based on the shared units file.

**Example (validation YAML):**
```yaml
field_validations:
  mass_unit:
    allowed_values_from: ../config/metric_units.yaml:mass
  mass_value:
    min_value: 0
    max_value: 1e12
```

---

## 5. Unit Conversion and Validation

- At runtime, all unit validation and conversion is performed using Pint (or a similar library), which is always in sync with the units file.
- The `metric_units.yaml` file is generated from Pint's registry, ensuring it is always up-to-date and standards-compliant.

---

## 6. Workflow
1. Generate `metric_units.yaml` from Pint and store it in a shared location (e.g., `config/metric_units.yaml`).
2. Entity and validation configs reference this file for allowed units.
3. All modules (extraction, validation, conversion, UI, etc.) use this file for static validation and Pint for runtime logic.

---

## 7. Handling Multiples (SI Prefixes) and Derived Units

- Multiples (e.g., kilo-, mega-, milli-, micro-, nano-) and derived units (e.g., acceleration, pressure) are included in `metric_units.yaml`.
- The units file can either explicitly list all common multiples and derived units or include a section for allowed prefixes and base unit compositions.
- Pint is used for all dynamic validation and conversion.

---

## 8. Benefits
- **Single source of truth** for all allowed units.
- **Consistency** across all modules and domains.
- **Easy to update** (regenerate from Pint if standards change).
- **Extensible** (add new units or prefixes centrally).
- **Industry-compliant** (follows SI/NIST standards).

---

## 9. Example File Locations
```
project-root/
  config/
    metric_units.yaml   # <--- Shared, generated from Pint
  common_dictionary/
    config/
      domains/
        chemistry/
          entity_config.yaml
          validation_config.yaml
  ... (other modules)
```

---

## 10. Example Workflow
1. Entity config defines `mass_value` and `mass_unit`.
2. Validation config says: `mass_unit` must be in `config/metric_units.yaml` under `mass`.
3. `metric_units.yaml` lists all allowed mass units (including multiples and derived units).
4. Conversion script loads `metric_units.yaml`, validates units, and performs conversions as needed using Pint.

---

## 11. References
- [SI Brochure (BIPM)](https://www.bipm.org/en/publications/si-brochure)
- [NIST Guide to SI Units](https://physics.nist.gov/cuu/Units/)
- [Pint Python Library](https://pint.readthedocs.io/en/stable/)

---

**End of Document** 