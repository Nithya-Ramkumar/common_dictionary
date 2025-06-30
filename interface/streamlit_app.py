import streamlit as st
import yaml
import json
import os
import tempfile
from pathlib import Path
import sys

# Add the src directory to the path
sys.path.append(str(Path(__file__).parent.parent / "src"))

from config.env_loader import EnvironmentLoader
from config_reconciliation import ConfigReconciliation

# Custom CSS for better styling
st.set_page_config(
    page_title="Common Dictionary - Config Reconciliation",
    page_icon="🔧",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better fonts and layout
st.markdown("""
<style>
    .main-header {
        font-size: 3rem !important;
        font-weight: bold !important;
        color: #1f77b4 !important;
        margin-bottom: 2rem !important;
    }
    .section-header {
        font-size: 2rem !important;
        font-weight: bold !important;
        color: #2c3e50 !important;
        margin-top: 2rem !important;
        margin-bottom: 1rem !important;
    }
    .error-box {
        background-color: #ffebee !important;
        border-left: 5px solid #f44336 !important;
        padding: 1rem !important;
        margin: 1rem 0 !important;
        border-radius: 5px !important;
    }
    .warning-box {
        background-color: #fff3e0 !important;
        border-left: 5px solid #ff9800 !important;
        padding: 1rem !important;
        margin: 1rem 0 !important;
        border-radius: 5px !important;
    }
    .success-box {
        background-color: #e8f5e8 !important;
        border-left: 5px solid #4caf50 !important;
        padding: 1rem !important;
        margin: 1rem 0 !important;
        border-radius: 5px !important;
    }
    .fix-suggestion {
        background-color: #f0f8ff !important;
        border: 2px solid #2196f3 !important;
        padding: 1.5rem !important;
        margin: 1rem 0 !important;
        border-radius: 8px !important;
    }
    .metric-card {
        background-color: #f8f9fa !important;
        padding: 1.5rem !important;
        border-radius: 10px !important;
        border: 1px solid #dee2e6 !important;
        text-align: center !important;
    }
    .metric-value {
        font-size: 2.5rem !important;
        font-weight: bold !important;
    }
    .metric-label {
        font-size: 1.2rem !important;
        color: #6c757d !important;
    }
    .config-file {
        background-color: #f8f9fa !important;
        border: 1px solid #dee2e6 !important;
        padding: 1rem !important;
        border-radius: 5px !important;
        margin: 0.5rem 0 !important;
    }
    .stButton > button {
        font-size: 1.2rem !important;
        padding: 0.75rem 1.5rem !important;
        border-radius: 8px !important;
    }
    .stSelectbox > div > div {
        font-size: 1.1rem !important;
    }
    .stRadio > div {
        font-size: 1.1rem !important;
    }
    .stExpander > div > div {
        font-size: 1.1rem !important;
    }
</style>
""", unsafe_allow_html=True)

def main():
    st.markdown('<h1 class="main-header">🔧 Configuration Reconciliation & Validation</h1>', unsafe_allow_html=True)
    st.markdown('<h2 class="section-header">Module 1: Entity Extraction and Enrichment - Configuration Validation</h2>', unsafe_allow_html=True)
    
    # Sidebar for configuration
    with st.sidebar:
        st.markdown('<h3 class="section-header">Configuration</h3>', unsafe_allow_html=True)
        
        # Environment selection
        environment = st.selectbox(
            "Environment",
            ["development", "production", "testing", "staging"],
            index=0
        )
        
        # Config source selection
        config_source = st.radio(
            "Configuration Source",
            ["Environment Default", "Upload Custom Files"],
            help="Use environment default configs or upload custom files for testing"
        )
        
        st.divider()
        
        # Environment info
        st.markdown('<h4 class="section-header">Environment Info</h4>', unsafe_allow_html=True)
        env_loader = EnvironmentLoader(environment=environment)
        st.write(f"**Environment:** {environment}")
        st.write(f"**Config Root:** {env_loader.get('CONFIG_ROOT', 'Not set')}")
        st.write(f"**Debug Mode:** {env_loader.get('DEBUG', False)}")
        
        # Show available config paths
        st.markdown('<h4 class="section-header">Config Paths</h4>', unsafe_allow_html=True)
        config_paths = {
            'Entity Config': env_loader.get('ENTITY_CONFIG'),
            'Source Mapping': env_loader.get('SOURCE_MAPPING'),
            'Conflict Resolution': env_loader.get('CONFLICT_RESOLUTION'),
            'Validation Config': env_loader.get('VALIDATION_CONFIG'),
            'Extraction Config': env_loader.get('EXTRACTION_CONFIG'),
        }
        
        for name, path in config_paths.items():
            if path and os.path.exists(path):
                st.success(f"✓ {name}")
            elif path:
                st.error(f"✗ {name}")
            else:
                st.warning(f"? {name}")
    
    # Main content area - Full width
    st.markdown('<h3 class="section-header">Configuration Validation</h3>', unsafe_allow_html=True)
    
    if config_source == "Environment Default":
        # Use environment default configs
        st.info("Using environment default configuration files")
        
        # Initialize reconciliation with environment
        reconciler = ConfigReconciliation(env_loader)
        
        # Run validation
        if st.button("🔍 Validate Configurations", type="primary", use_container_width=True):
            with st.spinner("Validating configurations..."):
                result = reconciler.validate_configs()
            
            display_validation_results(result)
            
    else:
        # Upload custom files
        st.info("Upload custom configuration files for testing")
        
        uploaded_files = {}
        file_labels = {
            'entity_config': 'Entity Configuration (YAML)',
            'source_mapping': 'Source Mapping (YAML)',
            'conflict_resolution': 'Conflict Resolution (YAML)',
            'validation_config': 'Validation Configuration (YAML)',
            'extraction_config': 'Extraction Configuration (YAML)'
        }
        
        # File upload section
        for file_key, label in file_labels.items():
            uploaded_file = st.file_uploader(
                label,
                type=['yaml', 'yml'],
                key=file_key,
                help=f"Upload {label.lower()}"
            )
            if uploaded_file:
                uploaded_files[file_key] = uploaded_file
        
        # Validation button
        if st.button("🔍 Validate Uploaded Configurations", type="primary", disabled=len(uploaded_files) == 0, use_container_width=True):
            if len(uploaded_files) > 0:
                with st.spinner("Validating uploaded configurations..."):
                    result = validate_uploaded_files(uploaded_files, env_loader)
                display_validation_results(result)
            else:
                st.warning("Please upload at least one configuration file")
    
    # Quick Actions section
    st.markdown('<h3 class="section-header">Quick Actions</h3>', unsafe_allow_html=True)
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        if st.button("🧪 Test with Incomplete Configs", use_container_width=True):
            st.info("This would test with intentionally broken configurations")
            # TODO: Implement test with incomplete configs
    
    with col2:
        if st.button("📤 Export Results", use_container_width=True):
            st.info("Export validation results to JSON")
            # TODO: Implement export functionality
    
    with col3:
        if st.button("👁️ View Current Configs", use_container_width=True):
            display_current_configs(env_loader)
    
    st.divider()
    
    # Help section
    st.markdown('<h3 class="section-header">Help</h3>', unsafe_allow_html=True)
    st.markdown("""
    **Configuration Files:**
    - **Entity Config:** Defines entity types and attributes
    - **Source Mapping:** Maps attributes to data sources
    - **Conflict Resolution:** Defines conflict resolution strategies
    - **Validation Config:** Defines validation rules
    - **Extraction Config:** Defines extraction parameters
    
    **Validation Checks:**
    - File existence and YAML syntax
    - Cross-reference consistency
    - Missing or orphaned mappings
    - Required field validation
    """)

def validate_uploaded_files(uploaded_files, env_loader):
    """Validate uploaded configuration files"""
    # Create temporary config paths
    config_paths = {}
    
    for file_key, uploaded_file in uploaded_files.items():
        # Save uploaded file to temporary location
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as tmp_file:
            content = uploaded_file.read().decode('utf-8')
            tmp_file.write(content)
            config_paths[file_key] = tmp_file.name
    
    # Fill missing paths with environment defaults
    default_paths = {
        'entity_config': env_loader.get('ENTITY_CONFIG'),
        'source_mapping': env_loader.get('SOURCE_MAPPING'),
        'conflict_resolution': env_loader.get('CONFLICT_RESOLUTION'),
        'validation_config': env_loader.get('VALIDATION_CONFIG'),
        'extraction_config': env_loader.get('EXTRACTION_CONFIG'),
    }
    
    for key, default_path in default_paths.items():
        if key not in config_paths and default_path:
            config_paths[key] = default_path
    
    # Initialize reconciliation with custom paths
    reconciler = ConfigReconciliation(env_loader, config_paths)
    return reconciler.validate_configs()

def display_validation_results(result):
    """Display validation results in a structured way with immediate fix suggestions"""
    st.markdown('<h3 class="section-header">Validation Results</h3>', unsafe_allow_html=True)
    
    # Overall status with big, prominent display
    if result['success']:
        st.markdown('<div class="success-box"><h2>✅ All configurations are valid and consistent!</h2></div>', unsafe_allow_html=True)
        st.balloons()
    else:
        st.markdown(f'<div class="error-box"><h2>❌ Validation failed with {len(result["errors"])} errors</h2></div>', unsafe_allow_html=True)
    
    # Environment info
    st.info(f"**Environment:** {result['environment']}")
    
    # Results tabs
    tab1, tab2, tab3, tab4 = st.tabs(["📊 Summary", "❌ Errors & Fixes", "⚠️ Warnings", "ℹ️ Information"])
    
    with tab1:
        # Summary statistics with custom styling
        st.markdown('<h4 class="section-header">Summary Statistics</h4>', unsafe_allow_html=True)
        
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.markdown(f'''
            <div class="metric-card">
                <div class="metric-value">{len(result['errors'])}</div>
                <div class="metric-label">Errors</div>
            </div>
            ''', unsafe_allow_html=True)
        with col2:
            st.markdown(f'''
            <div class="metric-card">
                <div class="metric-value">{len(result['warnings'])}</div>
                <div class="metric-label">Warnings</div>
            </div>
            ''', unsafe_allow_html=True)
        with col3:
            st.markdown(f'''
            <div class="metric-card">
                <div class="metric-value">{len(result['info'])}</div>
                <div class="metric-label">Info</div>
            </div>
            ''', unsafe_allow_html=True)
        with col4:
            status = "✅ Pass" if result['success'] else "❌ Fail"
            status_color = "#4caf50" if result['success'] else "#f44336"
            st.markdown(f'''
            <div class="metric-card">
                <div class="metric-value" style="color: {status_color};">{status}</div>
                <div class="metric-label">Status</div>
            </div>
            ''', unsafe_allow_html=True)
        
        # Config paths
        st.markdown('<h4 class="section-header">Configuration Paths</h4>', unsafe_allow_html=True)
        for config_type, path in result['config_paths'].items():
            if path and os.path.exists(path):
                st.markdown(f'<div class="config-file">✅ <strong>{config_type}:</strong> {path}</div>', unsafe_allow_html=True)
            elif path:
                st.markdown(f'<div class="config-file" style="border-color: #f44336;">❌ <strong>{config_type}:</strong> {path}</div>', unsafe_allow_html=True)
            else:
                st.markdown(f'<div class="config-file" style="border-color: #ff9800;">⚠️ <strong>{config_type}:</strong> Not set</div>', unsafe_allow_html=True)
    
    with tab2:
        if result['errors']:
            st.markdown('<h4 class="section-header">Errors Found</h4>', unsafe_allow_html=True)
            st.markdown("**Click on any error to see detailed information and suggested fixes.**")
            
            for i, error in enumerate(result['errors'], 1):
                # Get fix suggestion immediately
                fix_suggestion = get_fix_suggestion(error)
                related_file = get_related_config_file(error, result['config_paths'])
                
                # Create expandable section with error and fix
                with st.expander(f"❌ Error {i}: {error[:80]}...", expanded=True):
                    # Error display
                    st.markdown(f'<div class="error-box"><strong>Error:</strong> {error}</div>', unsafe_allow_html=True)
                    
                    # Immediate fix suggestion
                    if fix_suggestion:
                        st.markdown('<div class="fix-suggestion"><strong>🔧 Suggested Fix:</strong></div>', unsafe_allow_html=True)
                        st.markdown(fix_suggestion)
                    
                    # Show related config file if applicable
                    if related_file and os.path.exists(related_file):
                        st.markdown('<div class="fix-suggestion"><strong>📄 Related Configuration File:</strong></div>', unsafe_allow_html=True)
                        try:
                            with open(related_file, 'r') as f:
                                content = f.read()
                            st.code(content, language='yaml')
                        except Exception as e:
                            st.error(f"Error reading file: {e}")
                    
                    # Add action buttons for common fixes
                    if 'missing mappings' in error.lower():
                        if st.button(f"🔧 Show Mapping Fix Examples for Error {i}", key=f"mapping_{i}"):
                            show_mapping_fix_examples()
                    elif 'orphaned entity type' in error.lower():
                        if st.button(f"🔧 Show Orphaned Fix Examples for Error {i}", key=f"orphaned_{i}"):
                            show_orphaned_fix_examples()
                    elif 'yaml syntax error' in error.lower():
                        if st.button(f"🔧 Show YAML Fix Examples for Error {i}", key=f"yaml_{i}"):
                            show_yaml_fix_examples()
        else:
            st.markdown('<div class="success-box"><h3>✅ No errors found!</h3></div>', unsafe_allow_html=True)
    
    with tab3:
        if result['warnings']:
            st.markdown('<h4 class="section-header">Warnings</h4>', unsafe_allow_html=True)
            for i, warning in enumerate(result['warnings'], 1):
                with st.expander(f"⚠️ Warning {i}: {warning[:80]}...", expanded=False):
                    st.markdown(f'<div class="warning-box"><strong>Warning:</strong> {warning}</div>', unsafe_allow_html=True)
                    
                    # Provide context-specific suggestions
                    suggestion = get_warning_suggestion(warning)
                    if suggestion:
                        st.markdown('<div class="fix-suggestion"><strong>💡 Suggestion:</strong></div>', unsafe_allow_html=True)
                        st.markdown(suggestion)
        else:
            st.markdown('<div class="success-box"><h3>✅ No warnings found!</h3></div>', unsafe_allow_html=True)
    
    with tab4:
        if result['info']:
            st.markdown('<h4 class="section-header">Information</h4>', unsafe_allow_html=True)
            for i, info in enumerate(result['info'], 1):
                st.markdown(f'<div class="success-box">{i}. {info}</div>', unsafe_allow_html=True)
        else:
            st.info("No additional information")

def get_fix_suggestion(error):
    """Get context-specific fix suggestions for errors"""
    error_lower = error.lower()
    
    if 'yaml syntax error' in error_lower:
        return """
        **🔧 YAML Syntax Error Fix:**
        
        **Common Issues:**
        1. **Unescaped backslashes** in regex patterns
        2. **Missing quotes** around strings with special characters
        3. **Incorrect indentation** (use spaces, not tabs)
        4. **Invalid escape sequences**
        
        **Example Fix:**
        ```yaml
        # ❌ Incorrect:
        pattern: "^[a-zA-Z0-9\s\-\.]+$"
        
        # ✅ Correct:
        pattern: "^[a-zA-Z0-9\\\\s\\\\-\\\\.]+$"
        ```
        
        **Quick Actions:**
        - Check the specific line mentioned in the error
        - Ensure all backslashes are properly escaped (double backslash)
        - Use quotes around strings containing special characters
        - Verify indentation uses spaces, not tabs
        """
    
    elif 'missing mappings' in error_lower:
        # Extract entity type and attributes from error message
        if 'entity' in error_lower and 'missing mappings for' in error_lower:
            return """
        **🔧 Missing Mappings Fix:**
        
        **What this means:** Attributes defined in your entity config are missing from source priority or conflict resolution configs.
        
        **How to fix:**
        1. **Add missing attributes to source_mapping.yaml:**
        ```yaml
        Compound:
          name:
            - source: pubchem
              priority: 1
          formula:
            - source: pubchem
              priority: 1
          molecular_mass_value:
            - source: pubchem
              priority: 1
            - source: reaxys
              priority: 2
        ```
        
        2. **Add missing attributes to conflict_resolution.yaml:**
        ```yaml
        Compound:
          name:
            strategy: "highest_priority"
            fallback: "first_source"
          formula:
            strategy: "highest_priority"
            fallback: "first_source"
          molecular_mass_value:
            strategy: "average"
            fallback: "highest_priority"
        ```
        
        **Quick Actions:**
        - Copy the missing attributes from the error message
        - Add them to both source_mapping.yaml and conflict_resolution.yaml
        - Ensure attribute names match exactly (case-sensitive)
        """
        else:
            return """
        **🔧 Missing Mappings Fix:**
        
        **What this means:** Some attributes are missing from your source priority or conflict resolution configurations.
        
        **How to fix:**
        1. Add the missing attributes to source_mapping.yaml
        2. Add the missing attributes to conflict_resolution.yaml
        3. Ensure attribute names match exactly (case-sensitive)
        
        **Example:**
        ```yaml
        # In source_mapping.yaml:
        EntityType:
          missing_attribute:
            - source: pubchem
              priority: 1
        ```
        """
    
    elif 'orphaned entity type' in error_lower:
        return """
        **🔧 Orphaned Entity Type Fix:**
        
        **What this means:** Entity types referenced in source priority or conflict resolution configs don't exist in the entity config.
        
        **How to fix:**
        
        1. **Add missing entity to entity_config.yaml:**
        ```yaml
        chemistry:
          types:
            - name: MissingEntityType
              description: "Description of the missing entity"
              subdomain: general
              attributes:
                - name: name
                  type: string
                  required: true
                - name: description
                  type: string
                  required: false
        ```
        
        2. **OR Remove orphaned references:**
        ```yaml
        # Remove from source_mapping.yaml:
        # MissingEntityType:  # Remove this entire section
        #   attribute_name:
        #     - source: some_source
        ```
        
        **Quick Actions:**
        - Check if the entity type should exist (add it to entity_config.yaml)
        - Or remove references to non-existent entity types from other configs
        - Check for typos in entity type names
        """
    
    elif 'file not found' in error_lower:
        return """
        **🔧 File Not Found Fix:**
        
        **What this means:** A configuration file cannot be found at the specified path.
        
        **How to fix:**
        1. **Check file path:** Verify the file exists in the specified location
        2. **Check permissions:** Ensure you have read access to the file
        3. **Check environment variables:** Verify CONFIG_ROOT and other path variables are set correctly
        4. **Create missing files:** If the file should exist, create it with proper YAML structure
        
        **Example file structure:**
        ```
        config/
        └── domains/
            └── chemistry/
                ├── entity_config.yaml
                ├── source_mapping.yaml
                ├── conflict_resolution.yaml
                ├── validation_config.yaml
                └── extraction_config.yaml
        ```
        
        **Quick Actions:**
        - Navigate to the directory and check if files exist
        - Verify environment variable paths in your .env file
        - Create missing files with proper YAML structure
        """
    
    elif 'no entities defined' in error_lower:
        return """
        **🔧 No Entities Defined Fix:**
        
        **What this means:** Your entity_config.yaml file doesn't contain any entity definitions.
        
        **How to fix:**
        
        **Add entities to entity_config.yaml:**
        ```yaml
        chemistry:
          version: "1.0.0"
          domain: chemistry
          types:
            - name: Compound
              description: "Chemical compounds and molecules"
              subdomain: general
              attributes:
                - name: name
                  type: string
                  required: true
                - name: formula
                  type: string
                  required: true
        ```
        
        **Quick Actions:**
        - Add at least one entity type to your entity_config.yaml
        - Ensure the entity has a name and at least one attribute
        - Follow the proper YAML structure with chemistry.types
        """
    
    return """
    **🔧 General Fix Suggestion:**
    
    **What this means:** An error was detected in your configuration files.
    
    **How to fix:**
    1. **Read the error message carefully** - it usually indicates the specific issue
    2. **Check the related configuration file** - the error will show which file needs attention
    3. **Verify YAML syntax** - ensure proper indentation and escaping
    4. **Check cross-references** - ensure entity types and attributes are consistent across files
    
    **Common issues:**
    - Missing or incorrect YAML syntax
    - Inconsistent entity type names across files
    - Missing required attributes or mappings
    - Incorrect file paths or permissions
    """

def get_warning_suggestion(warning):
    """Get context-specific suggestions for warnings"""
    warning_lower = warning.lower()
    
    if 'no attributes' in warning_lower:
        return """
        **Suggestion:** Consider adding attributes to this entity type 
        to make it more useful for data extraction and validation.
        """
    
    elif 'no validation rules' in warning_lower:
        return """
        **Suggestion:** Add validation rules for required attributes 
        to ensure data quality and consistency.
        """
    
    return "No specific suggestion available for this warning."

def get_related_config_file(error, config_paths):
    """Get the related config file for an error"""
    error_lower = error.lower()
    
    if 'entity_config' in error_lower or 'entities' in error_lower:
        return config_paths.get('entity_config')
    elif 'source_mapping' in error_lower or 'source priority' in error_lower:
        return config_paths.get('source_mapping')
    elif 'conflict_resolution' in error_lower:
        return config_paths.get('conflict_resolution')
    elif 'validation_config' in error_lower:
        return config_paths.get('validation_config')
    elif 'extraction_config' in error_lower:
        return config_paths.get('extraction_config')
    
    return None

def show_yaml_fix_examples():
    """Show examples of common YAML fixes"""
    st.markdown("""
    ### Common YAML Syntax Fixes
    
    **1. Escaping Backslashes in Regex Patterns:**
    ```yaml
    # ❌ Incorrect:
    pattern: "^[a-zA-Z0-9\\s\\-\\.]+$"
    
    # ✅ Correct:
    pattern: "^[a-zA-Z0-9\\\\s\\\\-\\\\.]+$"
    ```
    
    **2. Proper String Quoting:**
    ```yaml
    # ❌ Incorrect:
    description: This is a description with: colon
    
    # ✅ Correct:
    description: "This is a description with: colon"
    ```
    
    **3. Correct Indentation:**
    ```yaml
    # ❌ Incorrect (tabs):
    entities:
    	- name: Entity
    		attributes:
    			- name: attr
    
    # ✅ Correct (spaces):
    entities:
      - name: Entity
        attributes:
          - name: attr
    ```
    """)

def show_mapping_fix_examples():
    """Show examples of mapping fixes"""
    st.markdown("""
    ### Missing Mappings Fix Examples
    
    **1. Add Missing Source Priority Mapping:**
    ```yaml
    # In source_mapping.yaml:
    Compound:
      molecular_mass_value:
        - source: pubchem
          priority: 1
        - source: reaxys
          priority: 2
      molecular_mass_unit:
        - source: pubchem
          priority: 1
    ```
    
    **2. Add Missing Conflict Resolution:**
    ```yaml
    # In conflict_resolution.yaml:
    Compound:
      molecular_mass_value:
        strategy: "highest_priority"
        fallback: "average"
      molecular_mass_unit:
        strategy: "most_common"
        fallback: "SI_default"
    ```
    """)

def show_orphaned_fix_examples():
    """Show examples of orphaned entity fixes"""
    st.markdown("""
    ### Orphaned Entity Type Fix Examples
    
    **1. Add Missing Entity to entity_config.yaml:**
    ```yaml
    # In entity_config.yaml:
    chemistry:
      types:
        - name: MissingEntityType
          description: "Description of the missing entity"
          subdomain: general
          attributes:
            - name: name
              type: string
              required: true
            - name: description
              type: string
              required: false
    ```
    
    **2. Remove Orphaned References:**
    ```yaml
    # Remove from source_mapping.yaml:
    # MissingEntityType:  # Remove this entire section
    #   attribute_name:
    #     - source: some_source
    ```
    """)

def display_current_configs(env_loader):
    """Display current configuration files"""
    st.subheader("Current Configuration Files")
    
    config_files = {
        'Entity Config': env_loader.get('ENTITY_CONFIG'),
        'Source Mapping': env_loader.get('SOURCE_MAPPING'),
        'Conflict Resolution': env_loader.get('CONFLICT_RESOLUTION'),
        'Validation Config': env_loader.get('VALIDATION_CONFIG'),
        'Extraction Config': env_loader.get('EXTRACTION_CONFIG'),
    }
    
    for name, path in config_files.items():
        if path and os.path.exists(path):
            with st.expander(f"📄 {name}"):
                try:
                    with open(path, 'r') as f:
                        content = f.read()
                    st.code(content, language='yaml')
                except Exception as e:
                    st.error(f"Error reading file: {e}")
        else:
            st.warning(f"📄 {name}: File not found or path not set")

if __name__ == "__main__":
    main() 