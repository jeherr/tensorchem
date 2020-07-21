import json
from tensorchem.molecules.molecule import MoleculeSet

meta_natoms = {}
opt_natoms = {}
try:
    with open("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/meta_smiles_natoms.txt", "r") as f:
        meta_natoms = json.loads(f.read())
    meta_natoms = {int(key): value for key, value in meta_natoms.items()}
except FileNotFoundError:
    with open("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/meta_msets_smiles.txt", "r") as f:
        lines = f.readlines()
        for line in lines:
            meta_mol = line.rstrip('\n')
            meta_mol_data = json.load(open(meta_mol, "r"))
            if len(meta_mol_data['atoms']) in meta_natoms.keys():
                meta_natoms[len(meta_mol_data['atoms'])].append(meta_mol)
            else:
                meta_natoms[len(meta_mol_data['atoms'])] = [meta_mol]
    with open("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/meta_smiles_natoms.txt", "w") as file:
        json.dump(meta_natoms, file)
try:
    with open("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/opt_smiles_natoms.txt", "r") as f:
        opt_natoms = json.loads(f.read())
    opt_natoms = {int(key): value for key, value in opt_natoms.items()}
except FileNotFoundError:
    with open ("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/opt_msets_smiles.txt", "r") as f:
        lines = f.readlines()
        for line in lines:
            opt_mol = line.rstrip('\n')
            opt_mol_data = json.load(open(opt_mol, "r"))
            if len(opt_mol_data['atoms']) in opt_natoms.keys():
                opt_natoms[len(opt_mol_data['atoms'])].append(opt_mol)
            else:
                opt_natoms[len(opt_mol_data['atoms'])] = [opt_mol]
    with open("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/opt_smiles_natoms.txt", "w") as file:
        json.dump(opt_natoms, file)

opt_matches = {}
meta_matches = {}
for n_atoms, meta_mols in meta_natoms.items():
    if n_atoms in opt_natoms.keys():
        opt_mols = opt_natoms[n_atoms]
        opt_msets = []
        meta_msets = []
        for opt_mol in opt_mols:
            opt_mset = MoleculeSet()
            opt_mset.load(opt_mol)
            opt_mset.filename = opt_mol
            opt_msets.append(opt_mset)
        for meta_mol in meta_mols:
            meta_mset = MoleculeSet()
            meta_mset.load(meta_mol)
            meta_mset.filename = meta_mol
            meta_msets.append(meta_mset)
        for meta_mset in meta_msets:
            matches = []
            for opt_mset in opt_msets:
                if meta_mset.identifiers['smiles'] == opt_mset.identifiers['smiles']:
                    matches.append(opt_mset.filename)
            meta_matches[meta_mset.filename] = matches
        for opt_mset in opt_msets:
            matches = []
            for meta_mset in meta_msets:
                if opt_mset.identifiers['smiles'] == meta_mset.identifiers['smiles']:
                    matches.append(meta_mset.filename)
            opt_matches[opt_mset.filename] = matches

with open("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/opt_smile_matches.txt", "w") as f:
    json.dump(opt_matches, f)
with open("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/meta_smile_matches.txt", "w") as f:
    json.dump(meta_matches, f)
