from typing import Dict, Any, List
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from config.env_loader import EnvironmentLoader
from .base_source import BaseSource
import logging
import datetime

logger = logging.getLogger("pubchem")

class PubChemSource(BaseSource):
    """
    PubChem API data source implementation for search-based extraction.
    Implements the BaseSource interface for the new extraction flow.
    """
    # Flexible attribute mapping with endpoint, format, and batchability
    pubchem_attr_map = {
        'name': {
            'pubchem_property': 'IUPACName',
            'endpoint': 'property',
            'format': 'json',
            'batchable': True
        },
        'smiles': {
            'pubchem_property': 'SMILES',
            'endpoint': 'property',
            'format': 'json',
            'batchable': True
        },
        'canonical_smiles': {
            'pubchem_property': 'CanonicalSMILES',
            'endpoint': 'property',
            'format': 'json',
            'batchable': True
        },
        'isomeric_smiles': {
            'pubchem_property': 'IsomericSMILES',
            'endpoint': 'property',
            'format': 'json',
            'batchable': True
        },
        'connectivity_smiles': {
            'pubchem_property': 'ConnectivitySMILES',
            'endpoint': 'property',
            'format': 'json',
            'batchable': True
        },
        'molecular_formula': {
            'pubchem_property': 'MolecularFormula',
            'endpoint': 'property',
            'format': 'json',
            'batchable': True
        },
        'molecular_weight': {
            'pubchem_property': 'MolecularWeight',
            'endpoint': 'property',
            'format': 'json',
            'batchable': True
        },
        'pubchem_cid': {
            'pubchem_property': 'CID',
            'endpoint': 'property',
            'format': 'json',
            'batchable': True
        },
        'inchi': {
            'pubchem_property': 'InChI',
            'endpoint': 'property',
            'format': 'json',
            'batchable': True
        },
        'inchi_key': {
            'pubchem_property': 'InChIKey',
            'endpoint': 'property',
            'format': 'json',
            'batchable': True
        },
        'synonyms': {
            'pubchem_property': None,
            'endpoint': 'synonyms',
            'format': 'json',
            'batchable': False
        },
        # Add more as needed
    }

    def __init__(self, config, env_loader=None, debug=False):
        super().__init__(config, env_loader, debug)
        self.connection = config['connection']
        self.session = requests.Session()
        self.base_url = self.connection['base_url']
        self.timeout = self.connection.get('timeout', 30)
        self.retries = self.connection.get('retries', 3)
        self.debug = debug

    def search(self, entity_type: str, filters: List[Dict[str, Any]], attributes: List[str], max_results: int) -> List[Dict[str, Any]]:
        """
        Perform search-based extraction from PubChem.
        Returns a list of dicts, each with all requested attributes and a _provenance key.
        """
        # Determine search term
        search_term = None
        for f in filters:
            if f.get('attribute') == 'name':
                search_term = f.get('values', [None])[0]
        if not search_term:
            search_term = entity_type
        cids = self.get_cids_by_category(search_term, max_results)
        logger.debug(f"[PubChem] Search term: {search_term}, CIDs: {cids}")
        if not cids:
            logger.warning(f"[PubChem] No CIDs found for search term '{search_term}'")
            return []
        # Group attributes by batchability
        batch_attrs = []
        single_attrs = []
        for attr in attributes:
            config = self.pubchem_attr_map.get(attr)
            if config:
                if config['batchable']:
                    batch_attrs.append((attr, config))
                else:
                    single_attrs.append((attr, config))
        results = []
        # Batch fetch for batchable attributes
        if batch_attrs:
            # Only request the 'SMILES' property for batch fetch
            prop_names = [cfg['pubchem_property'] for attr, cfg in batch_attrs if cfg['pubchem_property'] and cfg['pubchem_property'] != 'CID']
            prop_list = ','.join(prop_names)
            url = f"{self.base_url}/compound/cid/{','.join(cids)}/property/{prop_list}/JSON"
            try:
                resp = self.session.get(url, timeout=self.timeout)
                resp.raise_for_status()
                data = resp.json()
                for row in data.get('PropertyTable', {}).get('Properties', []):
                    result = {}
                    for attr, cfg in batch_attrs:
                        pubchem_field = cfg['pubchem_property']
                        # Do not try to extract 'CID' as a property, handle it separately
                        if pubchem_field and pubchem_field != 'CID':
                            value = row.get(pubchem_field)
                            result[attr] = value if value is not None else 'unavailable'
                    # Always include pubchem_cid if available
                    if 'pubchem_cid' in [attr for attr, _ in batch_attrs]:
                        result['pubchem_cid'] = row.get('CID', 'unavailable')
                    # Add provenance
                    result['_provenance'] = {
                        'source': 'pubchem',
                        'endpoint': 'property',
                        'search_term': search_term,
                        'cids': [row.get('CID', 'unavailable')],
                        'property_list': prop_names,
                        'url': url,
                        'timestamp': datetime.datetime.utcnow().isoformat() + 'Z'
                    }
                    results.append(result)
            except Exception as e:
                logger.error(f"[PubChem] Batch property fetch failed: {e}")
        # Single fetch for non-batchable attributes (e.g., synonyms)
        for attr, cfg in single_attrs:
            for cid in cids:
                url = f"{self.base_url}/compound/cid/{cid}/{cfg['endpoint']}/{cfg['format']}"
                try:
                    resp = self.session.get(url, timeout=self.timeout)
                    resp.raise_for_status()
                    value = None
                    if cfg['endpoint'] == 'synonyms' and cfg['format'] == 'json':
                        data = resp.json()
                        value = data.get('InformationList', {}).get('Information', [{}])[0].get('Synonym', [])
                    # Add more endpoint/format handling as needed
                    result = {attr: value if value is not None else 'unavailable'}
                    # Add provenance
                    result['_provenance'] = {
                        'source': 'pubchem',
                        'endpoint': cfg['endpoint'],
                        'search_term': search_term,
                        'cids': [cid],
                        'property_list': [cfg['pubchem_property']],
                        'url': url,
                        'timestamp': datetime.datetime.utcnow().isoformat() + 'Z'
                    }
                    results.append(result)
                except Exception as e:
                    logger.error(f"[PubChem] Single property fetch failed for {attr}, CID {cid}: {e}")
        return results

    def extract_by_key(self, entity_type: str, key: Any, attributes: List[str]) -> Dict[str, Any]:
        """
        PubChem does not support key-based extraction for this use case.
        Returns 'unavailable' for all requested attributes, with provenance.
        """
        result = {attr: 'unavailable' for attr in attributes}
        result['_provenance'] = {
            'source': 'pubchem',
            'endpoint': None,
            'search_term': None,
            'cids': [],
            'property_list': [],
            'url': None,
            'timestamp': datetime.datetime.utcnow().isoformat() + 'Z'
        }
        return result

    def get_cids_by_category(self, term: str, count: int = 20) -> List[str]:
        """
        Fetch CIDs from PubChem by category search term using ESearch.
        """
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        params = {
            'db': 'pccompound',
            'term': term,
            'retmax': count,
            'retmode': 'json'
        }
        try:
            resp = self.session.get(url, params=params, timeout=self.timeout)
            resp.raise_for_status()
            data = resp.json()
            cids = data.get('esearchresult', {}).get('idlist', [])
            logger.debug(f"[PubChem] ESearch URL: {url}, term: {term}, CIDs: {cids}")
            return cids
        except Exception as e:
            logger.error(f"[PubChem] ESearch failed for term '{term}': {e}")
            return [] 