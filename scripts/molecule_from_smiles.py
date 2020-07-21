"""
Creates a Molecule object from a SMILES string using RDKit
"""

from rdkit import Chem

from torch import ByteTensor
from tensorchem.molecules import Molecule


def from_smiles(smiles: str) -> 'Molecule':
    rdmol = Chem.MolFromSmiles(smiles)
    rdmol = Chem.AddHs(rdmol)
    new_mol = Molecule()
    new_mol.atoms = ByteTensor([atom.GetAtomicNum() for atom in rdmol.GetAtoms()])
    new_mol.identifiers = {'smiles': smiles}
    return new_mol


if __name__ == "__main__":
    smile = "c"
    mol = from_smiles(smile)
    print(mol.atoms)
