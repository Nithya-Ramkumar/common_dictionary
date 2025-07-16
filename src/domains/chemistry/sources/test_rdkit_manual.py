from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem


def test_rdkit(smiles):
    print(f"Testing SMILES: {smiles}")
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        print("  Could not parse SMILES.")
        return
    try:
        mw = Descriptors.MolWt(mol)
        formula = rdMolDescriptors.CalcMolFormula(mol)
        rot_bonds = Descriptors.NumRotatableBonds(mol)
        # For rings, use mol.GetRingInfo().NumRings() for compatibility
        rings = mol.GetRingInfo().NumRings()
        logp = Descriptors.MolLogP(mol)
        print(f"  Molecular Weight: {mw}")
        print(f"  Molecular Formula: {formula}")
        print(f"  Num Rotatable Bonds: {rot_bonds}")
        print(f"  Num Rings: {rings}")
        print(f"  LogP: {logp}")
        # 3D structure
        mol_with_h = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol_with_h, randomSeed=0xf00d)
        AllChem.UFFOptimizeMolecule(mol_with_h)
        conf = mol_with_h.GetConformer()
        coords = [list(conf.GetAtomPosition(i)) for i in range(mol_with_h.GetNumAtoms())]
        print(f"  3D Coordinates (first 3 atoms): {coords[:3]} ...")
    except Exception as e:
        print(f"  Error during RDKit computation: {e}")

if __name__ == "__main__":
    # Example: Replace with your actual SMILES or add more to test
    smiles_list = [
        "C(COCCOCCN)N.N",
        # Add more SMILES strings here for testing
    ]
    for smi in smiles_list:
        test_rdkit(smi)
        print("-" * 40) 