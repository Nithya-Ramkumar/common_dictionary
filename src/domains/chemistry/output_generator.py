import os
import json
import csv
import yaml
from typing import List, Dict, Any
import datetime
import logging

logger = logging.getLogger("output_generator")

class OutputGenerator:
    """
    Generates output files (JSON, CSV, TSV, HTML) for extracted entities using a schema.
    The schema is loaded from output_entity_schema.yaml and determines which fields to include and their order.
    A summary is written at the top of each output file.
    
    Output formats:
    - JSON/CSV/TSV: Flat attribute records (one row per attribute extraction)
    - HTML: CID-grouped display with per-attribute metadata
    """
    def __init__(self, schema_path: str):
        self.schema = self._load_schema(schema_path)

    def _load_schema(self, schema_path: str) -> Dict[str, List[str]]:
        with open(schema_path, 'r') as f:
            data = yaml.safe_load(f)
        return data.get('output_schemas', {})

    def generate_outputs(self, entities: List[Dict[str, Any]], summary: Dict[str, Any], output_dir: str, schema_name: str = 'default', formats: List[str] = ['json', 'csv', 'tsv']):
        os.makedirs(output_dir, exist_ok=True)
        fields = self.schema.get(schema_name, [])
        
        # For JSON/CSV/TSV: Use flat attribute records (one row per attribute extraction)
        # For HTML: Group by CID for visual display
        if 'json' in formats:
            self._write_json(entities, summary, output_dir)
        if 'csv' in formats:
            self._write_table(entities, summary, output_dir, fields, sep=',', ext='csv')
        if 'tsv' in formats:
            self._write_table(entities, summary, output_dir, fields, sep='\t', ext='tsv')
        if 'html' in formats:
            self._write_html(entities, summary, output_dir, fields)

    def _write_json(self, entities: List[Dict[str, Any]], summary: Dict[str, Any], output_dir: str):
        """Write JSON output as flat attribute records with per-attribute metadata."""
        out_path = os.path.join(output_dir, 'entities.json')
        output = {
            'summary': summary,
            'entities': entities  # Flat list of attribute records
        }
        with open(out_path, 'w') as f:
            json.dump(output, f, indent=2)
        logger.info(f"Wrote JSON to {out_path}")

    def _write_table(self, entities: List[Dict[str, Any]], summary: Dict[str, Any], output_dir: str, fields: List[str], sep: str, ext: str):
        """Write CSV/TSV output as flat attribute records (one row per attribute extraction)."""
        out_path = os.path.join(output_dir, f'entities.{ext}')
        with open(out_path, 'w', newline='') as f:
            # Write summary as comments
            for k, v in summary.items():
                f.write(f"# {k}: {v}\n")
            
            # Define columns for flat attribute records, including expanded provenance fields
            columns = [
                'cid', 'attribute', 'value', 'timestamp', 'confidence', 'source',
                'provenance_source', 'provenance_method', 'provenance_endpoint', 'provenance_input_smiles',
                'provenance_attributes', 'provenance_cids', 'provenance_property_list', 'provenance_url'
            ]
            writer = csv.DictWriter(f, fieldnames=columns, delimiter=sep)
            writer.writeheader()
            
            # Write one row per attribute extraction
            for entity in entities:
                provenance = entity.get('provenance', {})
                row = {
                    'cid': entity.get('cid', ''),
                    'attribute': entity.get('attribute', ''),
                    'value': entity.get('value', ''),
                    'timestamp': entity.get('timestamp', ''),
                    'confidence': entity.get('confidence', ''),
                    'source': entity.get('source', ''),
                    'provenance_source': provenance.get('source', ''),
                    'provenance_method': provenance.get('method', ''),
                    'provenance_endpoint': provenance.get('endpoint', ''),
                    'provenance_input_smiles': provenance.get('input_smiles', ''),
                    'provenance_attributes': provenance.get('attributes', ''),
                    'provenance_cids': provenance.get('cids', ''),
                    'provenance_property_list': provenance.get('property_list', ''),
                    'provenance_url': provenance.get('url', ''),
                }
                writer.writerow(row)
        logger.info(f"Wrote {ext.upper()} to {out_path}")

    def _format_provenance_for_csv(self, provenance: Dict[str, Any]) -> str:
        """Format provenance dict as a string for CSV output."""
        if not provenance:
            return ''
        
        source = provenance.get('source', 'Unknown')
        endpoint = provenance.get('endpoint', '')
        query = provenance.get('query', {})
        
        parts = [source]
        if endpoint:
            parts.append(f"endpoint={endpoint}")
        if query:
            query_parts = [f"{k}={v}" for k, v in query.items()]
            parts.append(f"query=({', '.join(query_parts)})")
        
        return ' | '.join(parts)

    def _write_html(self, entities: List[Dict[str, Any]], summary: Dict[str, Any], output_dir: str, fields: List[str]):
        """Write HTML output with CID-grouped display and per-attribute metadata."""
        out_path = os.path.join(output_dir, 'entities.html')
        timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        polymers_requested = summary.get('polymers_requested', 'N/A')
        polymers_extracted = summary.get('polymers_extracted', len(entities))
        errors = summary.get('field_errors', {})
        
        # Group entities by CID for visual grouping
        entities_by_cid = {}
        for entity in entities:
            cid = entity.get('cid')
            if not cid:
                continue
            if cid not in entities_by_cid:
                entities_by_cid[cid] = []
            entities_by_cid[cid].append(entity)
        
        html = [
            '<!DOCTYPE html>',
            '<html>',
            '<head>',
            '<meta charset="UTF-8">',
            '<title>Polymer Extraction Report</title>',
            '<style>',
            'body { font-family: Arial, sans-serif; margin: 2em; background: #fafbfc; }',
            'h1 { color: #2c3e50; }',
            'h2 { color: #34495e; margin-top: 2em; }',
            'h3 { color: #3498db; background: #ecf0f1; padding: 10px; border-radius: 5px; margin: 1em 0 0.5em 0; }',
            'table { border-collapse: collapse; width: 100%; margin-top: 1em; }',
            'th, td { border: 1px solid #ccc; padding: 8px; text-align: left; }',
            'th { background: #f0f0f0; }',
            '.missing { background: #ffe0e0; }',
            '.provenance { font-size: 0.9em; color: #666; word-break: break-all; white-space: pre-line; }',
            '.summary-list { list-style: none; padding: 0; }',
            '.summary-list li { margin-bottom: 0.3em; }',
            '.cid-section { margin-bottom: 2em; border: 1px solid #ddd; border-radius: 8px; padding: 1em; background: white; }',
            '.cid-header { background: #3498db; color: white; padding: 10px; margin: -1em -1em 1em -1em; border-radius: 8px 8px 0 0; }',
            '.attribute-table { margin-top: 1em; }',
            '.attribute-table th { background: #ecf0f1; font-weight: bold; }',
            '.attribute-table td { background: #f8f9fa; word-break: break-all; white-space: pre-line; }',
            '.nav-tabs { display: flex; border-bottom: 1px solid #ddd; margin-bottom: 1em; }',
            '.nav-tab { padding: 10px 20px; cursor: pointer; border: 1px solid #ddd; border-bottom: none; background: #f8f9fa; margin-right: 5px; border-radius: 5px 5px 0 0; }',
            '.nav-tab.active { background: white; border-bottom: 1px solid white; margin-bottom: -1px; }',
            '.tab-content { display: none; }',
            '.tab-content.active { display: block; }',
            '</style>',
            '<script>',
            'function showTab(tabName) {',
            '  // Hide all tab contents',
            '  var tabContents = document.getElementsByClassName("tab-content");',
            '  for (var i = 0; i < tabContents.length; i++) {',
            '    tabContents[i].classList.remove("active");',
            '  }',
            '  // Remove active class from all tabs',
            '  var tabs = document.getElementsByClassName("nav-tab");',
            '  for (var i = 0; i < tabs.length; i++) {',
            '    tabs[i].classList.remove("active");',
            '  }',
            '  // Show selected tab content and activate tab',
            '  document.getElementById(tabName).classList.add("active");',
            '  event.target.classList.add("active");',
            '}',
            '</script>',
            '</head>',
            '<body>',
            '<h1>Polymer Extraction Report</h1>',
            f'<div><b>Run timestamp:</b> {timestamp}</div>',
            '<h2>Summary</h2>',
            '<ul class="summary-list">',
            f'<li><b>Polymers requested:</b> {polymers_requested}</li>',
            f'<li><b>Polymers extracted:</b> {polymers_extracted}</li>',
            f'<li><b>Errors:</b> {sum(len(v) for v in errors.values()) if errors else 0}</li>',
            '</ul>',
            '<h2>Extracted Polymers (Grouped by CID)</h2>'
        ]
        
        # Add navigation tabs
        if entities_by_cid:
            html.append('<div class="nav-tabs">')
            cid_list = list(entities_by_cid.keys())
            for i, cid in enumerate(cid_list):
                active_class = ' active' if i == 0 else ''
                html.append(f'<div class="nav-tab{active_class}" onclick="showTab(\'cid-{cid}\')">CID {cid}</div>')
            html.append('</div>')
        
        # Create content for each CID
        for i, (cid, cid_entities) in enumerate(entities_by_cid.items()):
            active_class = ' active' if i == 0 else ''
            html.append(f'<div id="cid-{cid}" class="tab-content{active_class}">')
            html.append(f'<div class="cid-section">')
            html.append(f'<div class="cid-header"><h3>CID: {cid}</h3></div>')
            
            # Create attribute table for this CID
            html.append('<table class="attribute-table">')
            html.append('<tr><th>Attribute</th><th>Value</th><th>Source</th><th>Timestamp</th><th>Confidence</th><th>Provenance</th></tr>')
            
            # Display attributes for this CID
            for entity in cid_entities:
                attr = entity.get('attribute', '')
                value = entity.get('value', '')
                provenance = entity.get('provenance', {})
                timestamp = entity.get('timestamp', 'N/A')
                confidence = entity.get('confidence', 'N/A')
                # Format provenance
                prov_lines = []
                for k, v in provenance.items():
                    prov_lines.append(f"<b>{k}</b>: {v}")
                prov_html = '<br>'.join(prov_lines) if prov_lines else 'N/A'
                # Format value
                if value in [None, '', [], {}]:
                    value = 'N/A'
                    cell_class = ' class="missing"'
                else:
                    cell_class = ''
                html.append(f'<tr>')
                html.append(f'<td><strong>{attr}</strong></td>')
                html.append(f'<td{cell_class}>{value}</td>')
                html.append(f'<td class="provenance">{provenance.get("source", "N/A")}</td>')
                html.append(f'<td>{timestamp}</td>')
                html.append(f'<td>{confidence}</td>')
                html.append(f'<td>{prov_html}</td>')
                html.append(f'</tr>')
            html.append('</table>')
            html.append('</div></div>')
        html.append('</body></html>')
        with open(out_path, 'w') as f:
            f.write('\n'.join(html))
        logger.info(f"Wrote HTML to {out_path}") 