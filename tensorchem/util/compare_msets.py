import glob
import os
import json
from tensorchem.dataset.molecule import MoleculeSet


meta_natoms = {}
try:
    with open("/home/adriscoll/tensorchem/data/chemspider_data/chno_meta_natoms.txt", "r") as f:
        meta_natoms = json.loads(f.read())
    meta_natoms = {int(key): value for key, value in meta_natoms.items()}
except FileNotFoundError:
    for meta_mol in glob.glob ("/home/adriscoll/tensorchem/data/chemspider_data/chno_metamd_mset/*.mset"):
        meta_mol = json.loads(meta_mol.read())
        if len(meta_mol['atoms']) in meta_natoms.keys():
            meta_natoms[len(meta_mol['atoms'])].append(meta_mol)
        else:
            meta_natoms[len(meta_mol['atoms'])] = [meta_mol]
    with open("/home/adriscoll/tensorchem/data/chemspider_data/chno_meta_natoms.txt", "w") as f:
        json.dump(meta_natoms, f)

for mol in glob.glob("/home/adriscoll/tensorchem/data/chemspider_data/chno_opt_mset/*.mset"):
    mset = MoleculeSet()
    mset.load(mol)
    meta_mset = MoleculeSet()

    if mset.n_atoms in meta_natoms.keys():
        for meta in meta_natoms[mset.n_atoms]:
            meta_mset.load(meta)
            if mset.compare_hash(meta_mset):
                print("match")
                with open("/home/adriscoll/tensorchem/data/chemspider_data/matching_msets.txt", "a+") as f:
                    f.write("atoms match:\t opt: " + mol + "\t meta: " + meta + "\n")
