from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem
import json
import logging
from domains.chemistry.sources.rdkit_source import (
    get_three_d_structure,
    get_molecular_formula,
    get_average_molecular_weight,
    get_exact_molecular_weight,
    get_num_rotatable_bonds,
    get_num_rings,
    get_logp,
    has_substruct_match,
    get_degree_distribution,
    get_branching_index,
    get_morgan_fingerprint,
    get_maccs_fingerprint
)
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s | %(levelname)s | %(name)s | %(message)s'
)
logger = logging.getLogger("test_rdkit_manual")

def test_rdkit(smiles):
    logger.info(f"Testing SMILES: {smiles}")
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        logger.warning("  Could not parse SMILES.")
        return
    # --- Original property checks ---
    try:
        mw = Descriptors.MolWt(mol)
        formula = rdMolDescriptors.CalcMolFormula(mol)
        rot_bonds = Descriptors.NumRotatableBonds(mol)
        rings = mol.GetRingInfo().NumRings()
        logp = Descriptors.MolLogP(mol)
        mol_with_h = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol_with_h, randomSeed=0xf00d)
        AllChem.UFFOptimizeMolecule(mol_with_h)
        conf = mol_with_h.GetConformer()
        coords = [list(conf.GetAtomPosition(i)) for i in range(mol_with_h.GetNumAtoms())]
        logger.info(f"  Molecular Weight: {mw}")
        logger.info(f"  Molecular Formula: {formula}")
        logger.info(f"  Num Rotatable Bonds: {rot_bonds}")
        logger.info(f"  Num Rings: {rings}")
        logger.info(f"  LogP: {logp}")
        logger.info(f"  3D Coordinates (first 3 atoms): {coords[:3]} ...")
    except Exception as e:
        logger.error(f"  Error during RDKit computation: {e}")

    # --- Serializability checks for all extraction functions ---
    try:
        results = {
            "three_d_structure": get_three_d_structure(mol),
            "molecular_formula": get_molecular_formula(mol),
            "average_molecular_weight": get_average_molecular_weight(mol),
            "exact_molecular_weight": get_exact_molecular_weight(mol),
            "num_rotatable_bonds": get_num_rotatable_bonds(mol),
            "num_rings": get_num_rings(mol),
            "logp": get_logp(mol),
            "sulfonic_group": has_substruct_match(mol, "[SX4](=O)(=O)(O)O"),
            "carboxylic_group": has_substruct_match(mol, "[CX3](=O)[OX2H1]"),
            "phosphonic_group": has_substruct_match(mol, "[PX4](=O)(O)(O)O"),
            "degree_distribution": get_degree_distribution(mol),
            "branching_index": get_branching_index(mol),
            "morgan_fp": get_morgan_fingerprint(mol),
            "maccs_fp": get_maccs_fingerprint(mol),
        }
        for name, value in results.items():
            logger.debug(f"  {name}: {type(value)} -> {str(value)[:100]}{'...' if len(str(value))>100 else ''}")
            try:
                json.dumps({name: value})
                logger.debug(f"    ✓ {name} is JSON serializable")
            except Exception as e:
                logger.warning(f"    ✗ {name} is NOT JSON serializable: {e}")
    except Exception as e:
        logger.error(f"  Error during RDKit extraction function tests: {e}")

if __name__ == "__main__":
    smiles_list = [
        "CCO",
        # Problematic/complex SMILES for testing:
        "C(C1[C@@H]2[C@@H](C([C@@H](O1)O[C@H]3[C@@H](C([C@@H](O[C@H]4[C@@H](C([C@@H](O[C@H]5[C@@H](C([C@@H](O[C@H]6[C@@H](C([C@@H](O[C@H]7[C@@H](C([C@H](O[C@H]8[C@@H](C([C@H](O[C@H]9[C@@H](C([C@H](O1)OC9CO)O)O)OC8CO)O)O)OC7CO)O)O)OC6CO)O)O)OC5COCC(COCC9[C@@H]1[C@@H](C([C@@H](O9)O[C@H]2[C@@H](C([C@@H](O[C@H]5[C@@H](C([C@@H](O[C@H]6[C@@H](C([C@@H](O[C@H]7[C@@H](C([C@@H](O[C@H]8[C@@H](C([C@H](O[C@H]9[C@@H](C([C@H](O1)OC9CO)O)O)OC8CO)O)O)OC7CO)O)O)OC6CO)O)O)OC5CO)O)O)O)O)O)O)O)OC4CO)O)O)OC3CO)O)O)O)O)O",
        "C(COCCOCCN)N.N",
        "[2H]C1=C(C(=C(C(=C1C2(CCCCC2)C3=C(C(=C(C(=C3[2H])[2H])OC(=O)C)[2H])[2H])[2H])[2H])OC)[2H].Cl",
        "CC(=O)OC1=CC=C(C=C1)C2(CCCCC2)C3=CC=C(C=C3)OC.Cl",
        "CC1=C(C=CC(=C1)C(C)(C)C2=CC(=C(C=C2)OC(=O)C)C)OC.Cl",
        # Add more SMILES strings here for testing
    ]
    for smi in smiles_list:
        test_rdkit(smi)
        logger.info("-" * 40) 