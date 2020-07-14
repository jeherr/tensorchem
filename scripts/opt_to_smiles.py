from tensorchem.dataset.molecule import MoleculeSet as MSet
from tensorchem.util.xyz2mol import geom_to_smiles
from multiprocessing import Pool
from functools import partial
import json
import glob

def mset_to_smiles(max_atoms, filename):
    mset = MSet()
    print(filename)
    try:
        mset.load(filename)
    except json.JSONDecodeError:
        print("json decode error")
        return None
    if mset.n_atoms > max_atoms:
        return None
    geom = mset.get_min_geom
    try:
        charges = [float(atom.labels['wB97X-D.6-311g**.charges']) for atom in geom.atoms]
    except KeyError:
        charges = [float(atom.labels['wb97x-d.6-311gss.mulliken_charges']) for atom in geom.atoms]
    charge = sum(charges)/len(charges)
    try:
        smiles = geom_to_smiles(geom, charge)
    except:
        return None
    mset.identifiers.update({"smiles": smiles})
    mset.save(filename)
    return filename

if __name__ == "__main__":
    with open("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/meta_smiles_natoms.txt", "r") as f:
        meta_contents = json.loads(f.read())
    meta_natoms = [int(key) for key in list(meta_contents.keys())]
    func = partial(mset_to_smiles, max(meta_natoms))
    with open("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/opt_msets_smiles.txt", "w") as f:
        with Pool(32) as p:
            for i in range(1, 10):
                print(i)
                for file in p.imap_unordered(func, glob.iglob("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/opt/"+str(i)+"/*.mset"), 100):
                    if file is not None:
                        f.write(file+"\n")
            for x in ['b', 'br', 'i', 'p', 'p_new', 'se', 'si']:
                for file in p.imap_unordered(func, glob.iglob("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/opt_cs40/"+str(x)+"/*.mset"), 100):
                    if file is not None:
                        f.write(file+"\n")
