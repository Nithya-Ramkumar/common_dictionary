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
        # Substance-specific attributes
        'sid': {
            'pubchem_property': 'SID',
            'endpoint': 'property',
            'format': 'json',
            'batchable': True
        },
        'cas': {
            'pubchem_property': 'CAS',
            'endpoint': 'property',
            'format': 'json',
            'batchable': True
        },
        'brand_names': {
            'pubchem_property': 'TradeName',
            'endpoint': 'property',
            'format': 'json',
            'batchable': True
        },
        'source': {
            'pubchem_property': 'Source',
            'endpoint': 'property',
            'format': 'json',
            'batchable': True
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

    def search(self, entity_type: str, filters: List[Dict[str, Any]], attributes: List[str], max_results: int, all_values: Dict[str, List[Any]] = None) -> List[Dict[str, Any]]:
        """
        Perform search-based extraction from PubChem.
        Returns a list of dicts, each with all requested attributes and a _provenance key.
        """
        # If all_values is provided, use batch processing
        if all_values and 'name' in all_values:
            return self._batch_search_by_names(all_values['name'], attributes, max_results)
        
        # Original single-search logic (for backward compatibility)
        search_term = None
        for f in filters:
            if f.get('attribute') == 'name':
                search_term = f.get('values', [None])[0]
        if not search_term:
            search_term = entity_type
        
        # For substances (polymers), search for SIDs instead of CIDs
        if entity_type == 'Substance':
            sids = self.get_sids_by_category(search_term, max_results)
            logger.debug(f"[PubChem] Search term: {search_term}, SIDs: {sids}")
            if not sids:
                logger.warning(f"[PubChem] No SIDs found for search term '{search_term}'")
                return []
            return self._fetch_properties_for_sids(sids, attributes, search_term)
        else:
            # For compounds, use existing CID logic
            cids = self.get_cids_by_category(search_term, max_results)
            logger.debug(f"[PubChem] Search term: {search_term}, CIDs: {cids}")
            if not cids:
                logger.warning(f"[PubChem] No CIDs found for search term '{search_term}'")
                return []
            return self._fetch_properties_for_cids(cids, attributes, search_term)

    def _batch_search_by_names(self, names: List[str], attributes: List[str], max_results: int) -> List[Dict[str, Any]]:
        """
        Batch search for multiple names, then fetch properties for all CIDs/SIDs at once.
        """
        # For substances, use SIDs; for compounds, use CIDs
        if 'sid' in attributes or 'cas' in attributes or 'brand_names' in attributes or 'source' in attributes:
            # Substance search
            all_sids = set()
            search_terms_used = []
            
                    # Get SIDs for each name
        for name in names:
            # Add rate limiting delay
            self._add_rate_limiting_delay()
            sids = self.get_sids_by_category(name, max_results)
            if sids:
                all_sids.update(sids)
                search_terms_used.append(name)
            
            logger.debug(f"[PubChem] Batch search for {len(names)} names, found {len(all_sids)} unique SIDs")
            
            if not all_sids:
                logger.warning(f"[PubChem] No SIDs found for any of the search terms: {names}")
                return []
            
            # Fetch properties for all SIDs in one batch
            return self._fetch_properties_for_sids(list(all_sids), attributes, f"batch_search_{len(search_terms_used)}_terms")
        else:
            # Compound search (existing logic)
            all_cids = set()
            search_terms_used = []
            
            # Get CIDs for each name
            for name in names:
                cids = self.get_cids_by_category(name, max_results)
                if cids:
                    all_cids.update(cids)
                    search_terms_used.append(name)
            
            logger.debug(f"[PubChem] Batch search for {len(names)} names, found {len(all_cids)} unique CIDs")
            
            if not all_cids:
                logger.warning(f"[PubChem] No CIDs found for any of the search terms: {names}")
                return []
            
            # Fetch properties for all CIDs in one batch
            return self._fetch_properties_for_cids(list(all_cids), attributes, f"batch_search_{len(search_terms_used)}_terms")

    def _fetch_properties_for_cids(self, cids: List[str], attributes: List[str], search_term: str) -> List[Dict[str, Any]]:
        """
        Fetch properties for a list of CIDs using batch API calls.
        """
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

    def _fetch_properties_for_sids(self, sids: List[str], attributes: List[str], search_term: str) -> List[Dict[str, Any]]:
        """
        Fetch properties for a list of SIDs using PUG-View API for substances.
        """
        results = []
        
        # For substances, we need to use PUG-View API which doesn't support batch operations
        # So we fetch each SID individually
        for sid in sids:
            # Add rate limiting delay
            self._add_rate_limiting_delay()
            result = {}
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/substance/{sid}/JSON"
            try:
                resp = self.session.get(url, timeout=self.timeout)
                resp.raise_for_status()
                data = resp.json()
                
                # Extract properties from the PUG-View response
                record = data.get('Record', {})
                result['sid'] = sid
                result['name'] = record.get('RecordTitle', 'unavailable')
                
                # Extract SMILES and other properties from the record
                for section in record.get('Section', []):
                    if section.get('TOCHeading') == 'Names and Identifiers':
                        for subsection in section.get('Section', []):
                            if subsection.get('TOCHeading') == 'Computed Descriptors':
                                for info in subsection.get('Information', []):
                                    value = info.get('Value', {})
                                    if 'StringWithMarkup' in value:
                                        for item in value['StringWithMarkup']:
                                            string_val = item.get('String', '')
                                            if 'SMILES' in string_val:
                                                result['smiles'] = string_val.split('SMILES: ')[-1].strip()
                                            elif 'InChI=' in string_val:
                                                result['inchi'] = string_val
                                            elif 'InChIKey=' in string_val:
                                                result['inchi_key'] = string_val
                                            elif 'Molecular Formula:' in string_val:
                                                result['molecular_formula'] = string_val.split('Molecular Formula: ')[-1].strip()
                                            elif 'Molecular Weight:' in string_val:
                                                result['molecular_weight'] = string_val.split('Molecular Weight: ')[-1].strip()
                            elif subsection.get('TOCHeading') == 'Other Identifiers':
                                for info in subsection.get('Information', []):
                                    value = info.get('Value', {})
                                    if 'StringWithMarkup' in value:
                                        for item in value['StringWithMarkup']:
                                            string_val = item.get('String', '')
                                            if 'CAS:' in string_val:
                                                result['cas'] = string_val.split('CAS: ')[-1].strip()
                                            elif 'Trade Name:' in string_val:
                                                result['brand_names'] = string_val.split('Trade Name: ')[-1].strip()
                                            elif 'Source:' in string_val:
                                                result['source'] = string_val.split('Source: ')[-1].strip()
                
                # Add provenance
                result['_provenance'] = {
                    'source': 'pubchem',
                    'endpoint': 'pug_view',
                    'search_term': search_term,
                    'sids': [sid],
                    'property_list': list(result.keys()),
                    'url': url,
                    'timestamp': datetime.datetime.utcnow().isoformat() + 'Z'
                }
                
                # Mark missing attributes as unavailable
                for attr in attributes:
                    if attr not in result:
                        result[attr] = 'unavailable'
                
                results.append(result)
                
            except Exception as e:
                logger.error(f"[PubChem] PUG-View fetch failed for SID {sid}: {e}")
                # Add a result with unavailable values
                result = {attr: 'unavailable' for attr in attributes}
                result['sid'] = sid
                result['_provenance'] = {
                    'source': 'pubchem',
                    'endpoint': 'pug_view',
                    'search_term': search_term,
                    'sids': [sid],
                    'property_list': [],
                    'url': url,
                    'timestamp': datetime.datetime.utcnow().isoformat() + 'Z'
                }
                results.append(result)
        
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

    def get_sids_by_category(self, term: str, count: int = 20) -> List[str]:
        """
        Fetch SIDs from PubChem by category search term using ESearch.
        """
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        params = {
            'db': 'pcsubstance',
            'term': term,
            'retmax': count,
            'retmode': 'json'
        }
        try:
            resp = self.session.get(url, params=params, timeout=self.timeout)
            resp.raise_for_status()
            data = resp.json()
            sids = data.get('esearchresult', {}).get('idlist', [])
            logger.debug(f"[PubChem] ESearch URL: {url}, term: {term}, SIDs: {sids}")
            return sids
        except Exception as e:
            logger.error(f"[PubChem] ESearch failed for term '{term}': {e}")
            return []
        
    def _add_rate_limiting_delay(self):
        """
        Add a small delay to prevent rate limiting.
        """
        import time
        time.sleep(0.5)  # 500ms delay between requests 