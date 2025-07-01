# ... existing code ...
# TODO: Display and allow editing of subclass_of, subtypes, tags, enrichment rule details
# TODO: Show runtime feedback fields (confidence, source, version) for each attribute in review
# TODO: Visualize ontology/taxonomy and relationships (e.g., as a tree or graph)
# TODO: Load and display ontology.yaml for navigation and validation
# ... existing code ...

import yaml
import streamlit as st

def show_ontology_summary(ontology_path):
    with open(ontology_path, 'r') as f:
        ontology = yaml.safe_load(f)
    entity_names = [e['name'] for e in ontology.get('entity_types', [])]
    subdomains = ontology.get('subdomains', [])
    subclasses = {e['name']: e.get('subclass_of') for e in ontology.get('entity_types', []) if 'subclass_of' in e}
    subtypes = {e['name']: e.get('subtypes', []) for e in ontology.get('entity_types', []) if 'subtypes' in e}
    tags = {e['name']: e.get('tags', []) for e in ontology.get('entity_types', []) if 'tags' in e}
    st.subheader('Ontology Summary')
    st.markdown(f"**Subdomains:** {', '.join(subdomains) if subdomains else 'None'}")
    st.markdown(f"**Entity Types:** {', '.join(entity_names)}")
    if subclasses:
        st.markdown("**Subclass Relationships:**")
        for child, parent in subclasses.items():
            st.markdown(f"- {child} subclass_of {parent}")
    if subtypes:
        st.markdown("**Subtypes:**")
        for parent, children in subtypes.items():
            st.markdown(f"- {parent} has subtypes: {', '.join(children)}")
    if tags:
        st.markdown("**Tags:**")
        for entity, taglist in tags.items():
            st.markdown(f"- {entity}: {', '.join(taglist)}")

def main():
    # ... existing code ...
    st.title('Common Dictionary Config Validation & Ontology Viewer')
    # ... existing code ...
    # Show ontology summary
    show_ontology_summary('config/domains/chemistry/ontology.yaml')
    # TODO: Show validation results, errors, and allow navigation/fix
    # TODO: Display and allow editing of subclass_of, subtypes, tags, enrichment rule details
    # TODO: Show runtime feedback fields (confidence, source, version) for each attribute in review
    # TODO: Visualize ontology/taxonomy and relationships (e.g., as a tree or graph)
    # TODO: Load and display ontology.yaml for navigation and validation
    # ... existing code ...

if __name__ == "__main__":
    main()
# ... existing code ... 