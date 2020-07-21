from tensorchem.dataset.molecule import MoleculeSet as MSet
from tensorchem.util.xyz2mol import geom_to_smiles
from multiprocessing import Pool
from functools import partial
import json
import glob


def get_min_geom(geometries):
    energies = [geom.labels['potential'][0] for geom in geometries]
    return geometries[energies.index(min(energies))]


def mset_to_smiles(filename):
    mset = MSet()
    try:
        mset.load(filename)
    except json.JSONDecodeError:
        print("json decode error")
        return None
    geom = get_min_geom(mset.geometries)
    charge = 0
    try:
        smiles = geom_to_smiles(geom, charge)
    except:
        return None
    mset.identifiers.update({"smiles": smiles})
    mset.save(filename)
    return filename


def get_meta_smiles():
    with open("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/meta_msets_smiles.txt", "w") as f:
        with Pool(32) as p:
            for file in p.imap_unordered(mset_to_smiles, glob.glob("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/meta/*.mset")+glob.glob("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/aimd/*.mset"), 100):
                if file is not None:
                    f.write(file+"\n")


def get_opt_smiles():
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