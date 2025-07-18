from rdkit import Chem
from domains.chemistry.sources.rdkit_source import (
    get_molecular_formula,
    get_average_molecular_weight,
    get_exact_molecular_weight,
    get_num_rotatable_bonds,
    get_num_rings,
    get_logp,
    get_morgan_fingerprint,
    get_maccs_fingerprint
)

def test_one_molecule():
    smiles = "CCO"  # Ethanol
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, "Failed to parse SMILES"

    print("Molecular Formula:", get_molecular_formula(mol))
    print("Average Mol Wt:", get_average_molecular_weight(mol))
    print("Exact Mol Wt:", get_exact_molecular_weight(mol))
    print("Rotatable Bonds:", get_num_rotatable_bonds(mol))
    print("Rings:", get_num_rings(mol))
    print("LogP:", get_logp(mol))
    print("Morgan FP (first 10 bits):", get_morgan_fingerprint(mol)[:10])
    print("MACCS FP (first 10 bits):", get_maccs_fingerprint(mol)[:10])

if __name__ == "__main__":
    test_one_molecule()
