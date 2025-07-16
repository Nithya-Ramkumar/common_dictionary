from typing import Dict, Any, List
from .base_source import BaseSource
from config.env_loader import EnvironmentLoader

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem
except ImportError:
    Chem = None
    Descriptors = None
    rdMolDescriptors = None
    AllChem = None

class RDKitSource(BaseSource):
    """
    RDKit local cheminformatics extractor for key-based extraction (using SMILES).
    Implements the BaseSource interface for the new extraction flow.
    """
    attribute_map = {
        'three_d_structure': 'three_d_structure',
        'molecular_weight': 'MolWt',
        'molecular_formula': 'MolFormula',
        'num_rotatable_bonds': 'NumRotatableBonds',
        'num_rings': 'NumRings',
        'logp': 'MolLogP',
    }

    def __init__(self, config, env_loader=None, debug=False):
        super().__init__(config, env_loader, debug)
        self.available = Chem is not None
        # Debug RDKit environment if enabled
        debug_rdkit = False
        if env_loader is not None:
            # Try both upper and lower case for compatibility
            debug_rdkit = str(env_loader.get('DEBUG_RDKIT', False)).lower() == 'true' or \
                          str(env_loader.get('debug_rdkit', False)).lower() == 'true'
        if debug_rdkit:
            import sys
            if Chem is not None:
                print(f"[RDKit DEBUG] RDKit version: {getattr(Chem, '__version__', 'unknown')}")
            else:
                print("[RDKit DEBUG] RDKit is not available in this environment.")
            print(f"[RDKit DEBUG] Python executable: {sys.executable}")

    def search(self, entity_type: str, filters: List[Dict[str, Any]], attributes: List[str], max_results: int) -> List[Dict[str, Any]]:
        """
        RDKit does not support search-based extraction in this context.
        Returns an empty list.
        """
        return []

    def extract_by_key(self, entity_type: str, key: Any, attributes: List[str]) -> Dict[str, Any]:
        """
        Compute attributes from a SMILES string using RDKit.
        Args:
            entity_type: Name of the entity type
            key: The SMILES string
            attributes: List of attributes to extract
        Returns:
            Dict with all requested attributes (missing as 'unavailable')
        """
        result = {}
        if not self.available or not key:
            result = {attr: 'unavailable' for attr in attributes}
            # Add provenance
            result['_provenance'] = {
                'source': 'rdkit',
                'method': 'local',
                'input_smiles': key,
                'attributes': attributes,
                'timestamp': __import__('datetime').datetime.utcnow().isoformat() + 'Z'
            }
            return result
        mol = Chem.MolFromSmiles(key)
        if not mol:
            result = {attr: 'unavailable' for attr in attributes}
            # Add provenance
            result['_provenance'] = {
                'source': 'rdkit',
                'method': 'local',
                'input_smiles': key,
                'attributes': attributes,
                'timestamp': __import__('datetime').datetime.utcnow().isoformat() + 'Z'
            }
            return result
        for attr in attributes:
            try:
                if attr == 'three_d_structure':
                    # Generate 3D coordinates (as a placeholder, return a string or None)
                    mol_with_h = Chem.AddHs(mol)
                    AllChem.EmbedMolecule(mol_with_h, randomSeed=0xf00d)
                    AllChem.UFFOptimizeMolecule(mol_with_h)
                    conf = mol_with_h.GetConformer()
                    coords = [list(conf.GetAtomPosition(i)) for i in range(mol_with_h.GetNumAtoms())]
                    result[attr] = str(coords)
                elif attr == 'molecular_weight':
                    result[attr] = Descriptors.MolWt(mol)
                elif attr == 'molecular_formula':
                    result[attr] = rdMolDescriptors.CalcMolFormula(mol)
                elif attr == 'num_rotatable_bonds':
                    result[attr] = Descriptors.NumRotatableBonds(mol)
                elif attr == 'num_rings':
                    result[attr] = rdmolops.GetSSSR(mol)
                elif attr == 'logp':
                    result[attr] = Descriptors.MolLogP(mol)
                else:
                    result[attr] = 'unavailable'
            except Exception:
                result[attr] = 'unavailable'
        # Add provenance
        result['_provenance'] = {
            'source': 'rdkit',
            'method': 'local',
            'input_smiles': key,
            'attributes': attributes,
            'timestamp': __import__('datetime').datetime.utcnow().isoformat() + 'Z'
        }
        return result 