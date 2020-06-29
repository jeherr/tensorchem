import glob
import os
import json
from tensorchem.dataset.molecule import MoleculeSet

meta_natoms = {}
aimd_natoms = {}
try:
    with open("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/expanded_meta_smiles_natoms.txt", "r") as f:
        meta_natoms = json.loads(f.read())
    meta_natoms = {int(key): value for key, value in meta_natoms.items()}
except FileNotFoundError:
    with open("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/meta_msets_smiles.txt", "r") as f:
        lines = f.readlines()
        for line in lines:
            if 'meta' in line:
                meta_mol = line.rstrip('\n')
                meta_mol_data = json.load(open(meta_mol, "r"))
                if len(meta_mol_data['atoms']) in meta_natoms.keys():
                    meta_natoms[len(meta_mol_data['atoms'])].append(meta_mol)
                else:
                    meta_natoms[len(meta_mol_data['atoms'])] = [meta_mol]
    with open("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/expanded_meta_smiles_natoms.txt", "w") as file:
        json.dump(meta_natoms, file)
try:
    with open("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/expanded_aimd_smiles_natoms.txt", "r") as f:
        aimd_natoms = json.loads(f.read())
    aimd_natoms = {int(key): value for key, value in aimd_natoms.items()}
except FileNotFoundError:
    with open ("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/meta_msets_smiles.txt", "r") as f:
        lines = f.readlines()
        for line in lines:
            if 'aimd' in line:
                aimd_mol = line.rstrip('\n')
                aimd_mol_data = json.load(open(aimd_mol, "r"))
                if len(aimd_mol_data['atoms']) in aimd_natoms.keys():
                    aimd_natoms[len(aimd_mol_data['atoms'])].append(aimd_mol)
                else:
                    aimd_natoms[len(aimd_mol_data['atoms'])] = [aimd_mol]
    with open("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/expanded_aimd_smiles_natoms.txt", "w") as file:
        json.dump(aimd_natoms, file)

aimd_matches = {}
meta_matches = {}
for n_atoms, meta_mols in meta_natoms.items():
    if n_atoms in aimd_natoms.keys():
        aimd_mols = aimd_natoms[n_atoms]
        aimd_msets = []
        meta_msets = []
        for aimd_mol in aimd_mols:
            aimd_mset = MoleculeSet()
            aimd_mset.load(aimd_mol)
            aimd_mset.filename = aimd_mol
            aimd_msets.append(aimd_mset)
        for meta_mol in meta_mols:
            meta_mset = MoleculeSet()
            meta_mset.load(meta_mol)
            meta_mset.filename = meta_mol
            meta_msets.append(meta_mset)
        for meta_mset in meta_msets:
            matches = []
            for aimd_mset in aimd_msets:
                if meta_mset.identifiers['smiles'] == aimd_mset.identifiers['smiles']:
                    matches.append(aimd_mset.filename)
            meta_matches[meta_mset.filename] = matches
        for aimd_mset in aimd_msets:
            matches = []
            for meta_mset in meta_msets:
                if aimd_mset.identifiers['smiles'] == meta_mset.identifiers['smiles']:
                    matches.append(meta_mset.filename)
            aimd_matches[aimd_mset.filename] = matches

with open("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/expanded_aimd_smile_matches.txt", "w") as f:
    json.dump(aimd_matches, f)
with open("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/expanded_meta_smile_matches.txt", "w") as f:
    json.dump(meta_matches, f)
