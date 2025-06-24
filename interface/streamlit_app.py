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
from config.config_reconciliation import ConfigReconciliation

def main():
    st.set_page_config(
        page_title="Common Dictionary - Config Reconciliation",
        page_icon="🔧",
        layout="wide",
        initial_sidebar_state="expanded"
    )
    
    st.title("🔧 Configuration Reconciliation & Validation")
    st.markdown("**Module 1: Entity Extraction and Enrichment - Configuration Validation**")
    
    # Sidebar for configuration
    with st.sidebar:
        st.header("Configuration")
        
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
        st.subheader("Environment Info")
        env_loader = EnvironmentLoader(environment=environment)
        st.write(f"**Environment:** {environment}")
        st.write(f"**Config Root:** {env_loader.get('CONFIG_ROOT', 'Not set')}")
        st.write(f"**Debug Mode:** {env_loader.get('DEBUG', False)}")
        
        # Show available config paths
        st.subheader("Config Paths")
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
    
    # Main content area
    col1, col2 = st.columns([2, 1])
    
    with col1:
        st.header("Configuration Validation")
        
        if config_source == "Environment Default":
            # Use environment default configs
            st.info("Using environment default configuration files")
            
            # Initialize reconciliation with environment
            reconciler = ConfigReconciliation(env_loader)
            
            # Run validation
            if st.button("🔍 Validate Configurations", type="primary"):
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
            if st.button("🔍 Validate Uploaded Configurations", type="primary", disabled=len(uploaded_files) == 0):
                if len(uploaded_files) > 0:
                    with st.spinner("Validating uploaded configurations..."):
                        result = validate_uploaded_files(uploaded_files, env_loader)
                    display_validation_results(result)
                else:
                    st.warning("Please upload at least one configuration file")
    
    with col2:
        st.header("Quick Actions")
        
        # Test with incomplete configs
        if st.button("🧪 Test with Incomplete Configs"):
            st.info("This would test with intentionally broken configurations")
            # TODO: Implement test with incomplete configs
        
        # Export results
        if st.button("📤 Export Results"):
            st.info("Export validation results to JSON")
            # TODO: Implement export functionality
        
        # View current configs
        if st.button("👁️ View Current Configs"):
            display_current_configs(env_loader)
        
        st.divider()
        
        # Help section
        st.subheader("Help")
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
    """Display validation results in a structured way"""
    st.subheader("Validation Results")
    
    # Overall status
    if result['success']:
        st.success("✅ All configurations are valid and consistent!")
    else:
        st.error(f"❌ Validation failed with {len(result['errors'])} errors")
    
    # Environment info
    st.info(f"**Environment:** {result['environment']}")
    
    # Results tabs
    tab1, tab2, tab3, tab4 = st.tabs(["📊 Summary", "❌ Errors", "⚠️ Warnings", "ℹ️ Information"])
    
    with tab1:
        # Summary statistics
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Errors", len(result['errors']))
        with col2:
            st.metric("Warnings", len(result['warnings']))
        with col3:
            st.metric("Info", len(result['info']))
        with col4:
            st.metric("Status", "✅ Pass" if result['success'] else "❌ Fail")
        
        # Config paths
        st.subheader("Configuration Paths")
        for config_type, path in result['config_paths'].items():
            if path and os.path.exists(path):
                st.success(f"✓ {config_type}: {path}")
            elif path:
                st.error(f"✗ {config_type}: {path}")
            else:
                st.warning(f"? {config_type}: Not set")
    
    with tab2:
        if result['errors']:
            for i, error in enumerate(result['errors'], 1):
                st.error(f"{i}. {error}")
        else:
            st.success("No errors found!")
    
    with tab3:
        if result['warnings']:
            for i, warning in enumerate(result['warnings'], 1):
                st.warning(f"{i}. {warning}")
        else:
            st.success("No warnings found!")
    
    with tab4:
        if result['info']:
            for i, info in enumerate(result['info'], 1):
                st.info(f"{i}. {info}")
        else:
            st.info("No additional information")

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