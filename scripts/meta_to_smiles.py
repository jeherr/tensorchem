from tensorchem.dataset.molecule import MoleculeSet as MSet
from tensorchem.util.xyz2mol import geom_to_smiles
from multiprocessing import Pool
from functools import partial
import json
import glob

def mset_to_smiles(filename):
    mset = MSet()
    try:
        mset.load(filename)
    except json.JSONDecodeError:
        print("json decode error")
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
    with open("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/meta_msets_smiles.txt", "w") as f:
        with Pool(32) as p:
            for file in p.imap_unordered(mset_to_smiles, glob.glob("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/meta/*.mset")+glob.glob("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/aimd/*.mset"), 100):
                if file is not None:
                    f.write(file+"\n")
