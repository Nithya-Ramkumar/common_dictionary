from typing import Dict, Any, List
from .base_source import BaseSource
from config.env_loader import EnvironmentLoader
import logging


logger = logging.getLogger("rdkit")

# Graceful RDKit import
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, AllChem, rdmolops
except ImportError:
    logger.error("RDKit could not be imported. All RDKit features will be unavailable.")
    Chem = Descriptors = AllChem = rdmolops = None

try:
    from rdkit.Chem import rdMolDescriptors
except ImportError:
    logger.error("Failed to import rdMolDescriptors")
    rdMolDescriptors = None

try:
    from rdkit.Chem import MACCSkeys
except ImportError:
    logger.error("Failed to import MACCSkeys")
    MACCSkeys = None

try:
    from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
except ImportError:
    logger.error("Failed to import GetMorganGenerator")
    GetMorganGenerator = None


from rdkit import DataStructs

# ... etc for other modules

SMARTS_PATTERNS = {
    "sulfonic_group": "[SX4](=O)(=O)(O)O",
    "carboxylic_group": "[CX3](=O)[OX2H1]",
    "phosphonic_group": "[PX4](=O)(O)(O)O"
}

def safe(func):
    """Decorator to catch and log RDKit-related errors."""
    def wrapper(mol, *args, **kwargs):
        try:
            if mol is None:
                raise ValueError("Molecule is None")
            return func(mol, *args, **kwargs)
        except Exception as e:
            logger.error(f"RDKit error in {func.__name__}: {e}")
            return 'unavailable'
    return wrapper

@safe
def has_substruct_match(mol, pattern): return bool(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)))

@safe
def get_degree_distribution(mol): return [atom.GetDegree() for atom in mol.GetAtoms()]

@safe
def get_branching_index(mol): return sum(1 for atom in mol.GetAtoms() if atom.GetDegree() > 2)

@safe
def get_morgan_fingerprint(mol, radius=2, nBits=2048):
    if GetMorganGenerator is None:
        raise ImportError("GetMorganGenerator is not available")
    fp = GetMorganGenerator(radius, nBits).GetFingerprint(mol)
    return fp.ToBitString()

@safe
def get_maccs_fingerprint(mol):
    if MACCSkeys is None:
        raise ImportError("MACCSkeys is not available")
    return MACCSkeys.GenMACCSKeys(mol).ToBitString()

@safe
def get_exact_molecular_weight(mol): return Descriptors.ExactMolWt(mol)

@safe
def get_average_molecular_weight(mol): return Descriptors.MolWt(mol)

@safe
def get_three_d_structure(mol):
    mol_with_h = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_with_h, randomSeed=0xf00d)
    AllChem.UFFOptimizeMolecule(mol_with_h)
    conf = mol_with_h.GetConformer()
    return [list(conf.GetAtomPosition(i)) for i in range(mol_with_h.GetNumAtoms())]

@safe
def get_molecular_formula(mol):
    if rdMolDescriptors is None:
        raise ImportError("rdMolDescriptors is not available")
    return rdMolDescriptors.CalcMolFormula(mol)

@safe
def get_num_rotatable_bonds(mol): return Descriptors.NumRotatableBonds(mol)

@safe
def get_num_rings(mol): return int(mol.GetRingInfo().NumRings())

@safe
def get_logp(mol): return Descriptors.MolLogP(mol)

RDKit_ATTRIBUTE_FUNCTIONS = {
    "three_d_structure": get_three_d_structure,
    "molecular_formula": get_molecular_formula,
    "average_molecular_weight": get_average_molecular_weight,
    "exact_molecular_weight": get_exact_molecular_weight,
    "num_rotatable_bonds": get_num_rotatable_bonds,
    "num_rings": get_num_rings,
    "logp": get_logp,
    "sulfonic_group": lambda mol: has_substruct_match(mol, SMARTS_PATTERNS["sulfonic_group"]),
    "carboxylic_group": lambda mol: has_substruct_match(mol, SMARTS_PATTERNS["carboxylic_group"]),
    "phosphonic_group": lambda mol: has_substruct_match(mol, SMARTS_PATTERNS["phosphonic_group"]),
    "degree_distribution": get_degree_distribution,
    "branching_index": get_branching_index,
    "morgan_fp": get_morgan_fingerprint,
    "maccs_fp": get_maccs_fingerprint
}

class RDKitSource(BaseSource):
    """
    Robust RDKit extractor for config-based and SMILES-driven property extraction.
    """

    def __init__(self, config, env_loader=None, debug=False):
        super().__init__(config, env_loader, debug)
        self.available = Chem is not None
        if debug or (env_loader and env_loader.get_debug_flag("rdkit")):
            logger.setLevel(logging.DEBUG)
            logger.debug(f"RDKit availability: {self.available}")
            if self.available:
                logger.debug(f"RDKit version: {getattr(Chem, '__version__', 'unknown')}")

    def search(self, entity_type: str, filters: List[Dict[str, Any]], attributes: List[str], max_results: int) -> List[Dict[str, Any]]:
        return []

    def extract_by_key(self, entity_type: str, key: Any, attributes: List[Any]) -> Dict[str, Any]:
        logger.debug(f"RDKit.extract_by_key: entity_type={entity_type}, key={key}, attributes={attributes}")
        result = {}

        if not self.available or not key:
            logger.error(f"RDKit unavailable or missing key: {key}")
            return self._fallback_result(attributes, key)

        mol = Chem.MolFromSmiles(key)
        if not mol:
            logger.error(f"RDKit failed to parse SMILES: {key}")
            return self._fallback_result(attributes, key)

        for attr in attributes:
            name = attr['name'] if isinstance(attr, dict) else attr
            method = attr.get('method') if isinstance(attr, dict) else None
            params = attr.get('params', {}) if isinstance(attr, dict) else {}
            pattern = attr.get('pattern') if isinstance(attr, dict) else None

            try:
                if method == 'smarts' and pattern:
                    logger.debug(f"Extracting SMARTS attribute '{name}' with pattern '{pattern}'")
                    result[name] = has_substruct_match(mol, pattern)
                elif method == 'fingerprint' and name in RDKit_ATTRIBUTE_FUNCTIONS:
                    logger.debug(f"Extracting fingerprint '{name}' with params {params}")
                    result[name] = RDKit_ATTRIBUTE_FUNCTIONS[name](mol, **params)
                elif name in RDKit_ATTRIBUTE_FUNCTIONS:
                    logger.debug(f"Extracting property '{name}'")
                    result[name] = RDKit_ATTRIBUTE_FUNCTIONS[name](mol)
                else:
                    logger.warning(f"Unknown RDKit attribute '{name}' for SMILES {key}")
                    result[name] = 'unavailable'
            except Exception as e:
                logger.error(f"RDKit exception for {name} and SMILES '{key}': {e}")
                result[name] = 'unavailable'

        result['_provenance'] = {
            'source': 'rdkit',
            'method': 'local',
            'input_smiles': key,
            'attributes': attributes,
            'timestamp': __import__('datetime').datetime.utcnow().isoformat() + 'Z'
        }
        return result

    def _fallback_result(self, attributes, key):
        result = {
            (attr['name'] if isinstance(attr, dict) else attr): 'unavailable'
            for attr in attributes
        }
        result['_provenance'] = {
            'source': 'rdkit',
            'method': 'local',
            'input_smiles': key,
            'attributes': attributes,
            'timestamp': __import__('datetime').datetime.utcnow().isoformat() + 'Z'
        }
        return result
