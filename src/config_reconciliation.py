import os
import yaml
import json
from collections import defaultdict
from typing import Dict, List, Any, Optional, Tuple
from pathlib import Path
from config.env_loader import EnvironmentLoader

class ConfigReconciliation:
    """Enhanced configuration reconciliation and validation system"""
    
    def __init__(self, env_loader: Optional[EnvironmentLoader] = None, config_paths: Optional[Dict[str, str]] = None):
        """
        Initialize config reconciliation
        
        Args:
            env_loader: EnvironmentLoader instance (uses default if None)
            config_paths: Optional custom config paths to override environment
        """
        self.env_loader = env_loader or EnvironmentLoader()
        self.config_paths = config_paths or self._get_default_config_paths()
        self.errors = []
        self.warnings = []
        self.info = []
        
    def _get_default_config_paths(self) -> Dict[str, str]:
        """Get default config paths from environment"""
        return {
            'entity_config': self.env_loader.get('ENTITY_CONFIG'),
            'source_mapping': self.env_loader.get('SOURCE_MAPPING'),
            'conflict_resolution': self.env_loader.get('CONFLICT_RESOLUTION'),
            'validation_config': self.env_loader.get('VALIDATION_CONFIG'),
            'extraction_config': self.env_loader.get('EXTRACTION_CONFIG'),
        }
    
    def load_yaml(self, path: str) -> Optional[Dict[str, Any]]:
        """Load and validate YAML file"""
        try:
            if not os.path.exists(path):
                self.errors.append(f"File not found: {path}")
                return None
            
            with open(path, 'r') as f:
                data = yaml.safe_load(f)
                if data is None:
                    self.errors.append(f"Empty or invalid YAML file: {path}")
                    return None
                return data
        except yaml.YAMLError as e:
            self.errors.append(f"YAML syntax error in {path}: {str(e)}")
            return None
        except Exception as e:
            self.errors.append(f"Error loading {path}: {str(e)}")
            return None
    
    def get_entity_attributes(self, entity_config: Dict[str, Any]) -> Dict[str, set]:
        """Extract entity attributes from entity config"""
        result = defaultdict(set)
        
        # Handle the actual config structure: chemistry.types
        if 'chemistry' in entity_config and 'types' in entity_config['chemistry']:
            entities = entity_config['chemistry']['types']
        # Fallback to the expected structure: entities
        elif 'entities' in entity_config:
            entities = entity_config['entities']
        else:
            self.warnings.append("No entities found in entity_config.yaml")
            return result
        
        if not entities:
            self.warnings.append("No entities found in entity_config.yaml")
            return result
            
        for entity in entities:
            entity_type = entity.get('name')
            if not entity_type:
                self.warnings.append("Entity without name found in entity_config.yaml")
                continue
                
            attributes = entity.get('attributes', [])
            for attr in attributes:
                attr_name = attr.get('name')
                if attr_name:
                    result[entity_type].add(attr_name)
                else:
                    self.warnings.append(f"Attribute without name found in entity {entity_type}")
        
        return result
    
    def get_source_priority_mappings(self, source_mapping: Dict[str, Any]) -> Dict[str, set]:
        """Extract source priority mappings"""
        result = defaultdict(set)
        
        for entity_type, mappings in source_mapping.items():
            if isinstance(mappings, dict):
                for attr_name, source_config in mappings.items():
                    result[entity_type].add(attr_name)
            elif isinstance(mappings, list):
                for mapping in mappings:
                    if isinstance(mapping, dict) and 'attribute' in mapping:
                        result[entity_type].add(mapping['attribute'])
        
        return result
    
    def get_conflict_resolution_mappings(self, conflict_resolution: Dict[str, Any]) -> Dict[str, set]:
        """Extract conflict resolution mappings"""
        result = defaultdict(set)
        
        for entity_type, strategies in conflict_resolution.items():
            if isinstance(strategies, dict):
                for attr_name, strategy in strategies.items():
                    result[entity_type].add(attr_name)
            elif isinstance(strategies, list):
                for strategy in strategies:
                    if isinstance(strategy, dict) and 'attribute' in strategy:
                        result[entity_type].add(strategy['attribute'])
        
        return result
    
    def validate_entity_config(self, entity_config: Dict[str, Any]) -> bool:
        """Validate entity config structure and content"""
        if not entity_config:
            self.errors.append("Entity config is empty")
            return False
        
        # Handle the actual config structure: chemistry.types
        if 'chemistry' in entity_config and 'types' in entity_config['chemistry']:
            entities = entity_config['chemistry']['types']
        # Fallback to the expected structure: entities
        elif 'entities' in entity_config:
            entities = entity_config['entities']
        else:
            self.errors.append("No entities defined in entity config")
            return False
            
        if not entities:
            self.errors.append("No entities defined in entity config")
            return False
            
        for i, entity in enumerate(entities):
            if not isinstance(entity, dict):
                self.errors.append(f"Entity {i} is not a dictionary")
                continue
                
            name = entity.get('name')
            if not name:
                self.errors.append(f"Entity {i} has no name")
                continue
                
            attributes = entity.get('attributes', [])
            if not attributes:
                self.warnings.append(f"Entity '{name}' has no attributes")
                continue
                
            for j, attr in enumerate(attributes):
                if not isinstance(attr, dict):
                    self.errors.append(f"Attribute {j} in entity '{name}' is not a dictionary")
                    continue
                    
                attr_name = attr.get('name')
                if not attr_name:
                    self.errors.append(f"Attribute {j} in entity '{name}' has no name")
                    continue
                    
                # Validate required fields
                if attr.get('required', False) and not attr.get('validation'):
                    self.warnings.append(f"Required attribute '{attr_name}' in entity '{name}' has no validation rules")
        
        return len([e for e in self.errors if 'entity config' in e.lower()]) == 0
    
    def validate_source_mapping(self, source_mapping: Dict[str, Any]) -> bool:
        """Validate source mapping structure"""
        if not source_mapping:
            self.errors.append("Source mapping is empty")
            return False
            
        for entity_type, mappings in source_mapping.items():
            if not mappings:
                self.warnings.append(f"Entity '{entity_type}' has no source mappings")
                continue
                
            if isinstance(mappings, dict):
                for attr_name, source_config in mappings.items():
                    if not source_config:
                        self.warnings.append(f"No source configuration for '{attr_name}' in entity '{entity_type}'")
                    elif isinstance(source_config, list):
                        if not source_config:
                            self.warnings.append(f"Empty source list for '{attr_name}' in entity '{entity_type}'")
        
        return len([e for e in self.errors if 'source mapping' in e.lower()]) == 0
    
    def validate_conflict_resolution(self, conflict_resolution: Dict[str, Any]) -> bool:
        """Validate conflict resolution structure"""
        if not conflict_resolution:
            self.errors.append("Conflict resolution is empty")
            return False
            
        for entity_type, strategies in conflict_resolution.items():
            if not strategies:
                self.warnings.append(f"Entity '{entity_type}' has no conflict resolution strategies")
                continue
                
            if isinstance(strategies, dict):
                for attr_name, strategy in strategies.items():
                    if not strategy:
                        self.warnings.append(f"No conflict resolution strategy for '{attr_name}' in entity '{entity_type}'")
        
        return len([e for e in self.errors if 'conflict resolution' in e.lower()]) == 0
    
    def validate_cross_references(self, entity_attrs: Dict[str, set], 
                                source_priority_attrs: Dict[str, set], 
                                conflict_resolution_attrs: Dict[str, set]) -> bool:
        """Validate cross-references between config files"""
        all_entity_types = set(entity_attrs.keys())
        
        # Check for missing mappings in source priority
        for entity_type in entity_attrs:
            missing_in_sp = entity_attrs[entity_type] - source_priority_attrs.get(entity_type, set())
            if missing_in_sp:
                self.errors.append(f"[Source Priority] Entity '{entity_type}': Missing mappings for {sorted(missing_in_sp)}")
            
            # Check for missing mappings in conflict resolution
            missing_in_cr = entity_attrs[entity_type] - conflict_resolution_attrs.get(entity_type, set())
            if missing_in_cr:
                self.errors.append(f"[Conflict Resolution] Entity '{entity_type}': Missing mappings for {sorted(missing_in_cr)}")
        
        # Check for orphaned mappings in source priority
        for entity_type in source_priority_attrs:
            if entity_type not in all_entity_types:
                self.errors.append(f"[Source Priority] Orphaned entity type: '{entity_type}'")
                continue
                
            orphaned = source_priority_attrs[entity_type] - entity_attrs.get(entity_type, set())
            if orphaned:
                self.errors.append(f"[Source Priority] Entity '{entity_type}': Orphaned mappings {sorted(orphaned)}")
        
        # Check for orphaned mappings in conflict resolution
        for entity_type in conflict_resolution_attrs:
            if entity_type not in all_entity_types:
                self.errors.append(f"[Conflict Resolution] Orphaned entity type: '{entity_type}'")
                continue
                
            orphaned = conflict_resolution_attrs[entity_type] - entity_attrs.get(entity_type, set())
            if orphaned:
                self.errors.append(f"[Conflict Resolution] Entity '{entity_type}': Orphaned mappings {sorted(orphaned)}")
        
        return len(self.errors) == 0
    
    def validate_configs(self) -> Dict[str, Any]:
        """Main validation method"""
        self.errors = []
        self.warnings = []
        self.info = []
        
        # Load all config files
        entity_config = self.load_yaml(self.config_paths['entity_config'])
        source_mapping = self.load_yaml(self.config_paths['source_mapping'])
        conflict_resolution = self.load_yaml(self.config_paths['conflict_resolution'])
        validation_config = self.load_yaml(self.config_paths['validation_config'])
        extraction_config = self.load_yaml(self.config_paths['extraction_config'])
        
        # Validate individual config files
        if entity_config:
            self.validate_entity_config(entity_config)
        if source_mapping:
            self.validate_source_mapping(source_mapping)
        if conflict_resolution:
            self.validate_conflict_resolution(conflict_resolution)
        
        # Validate cross-references if all required files are loaded
        if entity_config and source_mapping and conflict_resolution:
            entity_attrs = self.get_entity_attributes(entity_config)
            source_priority_attrs = self.get_source_priority_mappings(source_mapping)
            conflict_resolution_attrs = self.get_conflict_resolution_mappings(conflict_resolution)
            
            self.validate_cross_references(entity_attrs, source_priority_attrs, conflict_resolution_attrs)
            
            # Add summary info
            self.info.append(f"Found {len(entity_attrs)} entity types")
            for entity_type, attrs in entity_attrs.items():
                self.info.append(f"  - {entity_type}: {len(attrs)} attributes")
        
        # Prepare result
        result = {
            'success': len(self.errors) == 0,
            'errors': self.errors,
            'warnings': self.warnings,
            'info': self.info,
            'config_paths': self.config_paths,
            'environment': self.env_loader.get('ENVIRONMENT', 'unknown')
        }
        
        return result
    
    def get_validation_summary(self) -> str:
        """Get a formatted validation summary"""
        result = self.validate_configs()
        
        summary = f"Configuration Validation Report\n"
        summary += f"Environment: {result['environment']}\n"
        summary += f"Config Paths: {result['config_paths']}\n\n"
        
        if result['info']:
            summary += "Information:\n"
            for info in result['info']:
                summary += f"  ✓ {info}\n"
            summary += "\n"
        
        if result['warnings']:
            summary += "Warnings:\n"
            for warning in result['warnings']:
                summary += f"  ⚠ {warning}\n"
            summary += "\n"
        
        if result['errors']:
            summary += "Errors:\n"
            for error in result['errors']:
                summary += f"  ✗ {error}\n"
            summary += "\n"
        
        if result['success']:
            summary += "✅ All configurations are valid and consistent!"
        else:
            summary += f"❌ Validation failed with {len(result['errors'])} errors"
        
        return summary

def main():
    """Command-line interface"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Validate configuration consistency')
    parser.add_argument('--env', default='development', help='Environment to use')
    parser.add_argument('--config-dir', help='Custom config directory')
    parser.add_argument('--output', help='Output results to JSON file')
    
    args = parser.parse_args()
    
    # Initialize with environment
    env_loader = EnvironmentLoader(environment=args.env)
    reconciler = ConfigReconciliation(env_loader)
    
    # Run validation
    result = reconciler.validate_configs()
    
    # Print summary
    print(reconciler.get_validation_summary())
    
    # Save to file if requested
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(result, f, indent=2)
        print(f"\nResults saved to {args.output}")

if __name__ == "__main__":
    main() 