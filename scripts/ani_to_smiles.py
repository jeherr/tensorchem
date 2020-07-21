from tensorchem.molecules.molecule import MoleculeSet as MSet
from tensorchem.util.xyz2mol import geom_to_smiles
from multiprocessing import Pool
import glob

def mset_to_smiles(filename):
    mset = MSet()
    mset.load(filename)
    geom = mset.get_min_geom
    try:
        cm5_charges = [float(atom.labels['wb97x_dz.cm5_charges']) for atom in geom.atoms]
        hirshfeld_charges = [float(atom.labels['wb97x_dz.hirshfeld_charges']) for atom in geom.atoms]
        charge = ((sum(cm5_charges)/len(cm5_charges)) + (sum(hirshfeld_charges)/len(hirshfeld_charges))) / 2
    except:
        charge = 0
    try:
        smiles = geom_to_smiles(geom, charge)
    except:
        return None
    mset.identifiers.update({"smiles": smiles})
    mset.save(filename)
    return filename

if __name__ == "__main__":
    with open("/mnt/sdb1/adriscoll/ani1x-data/ani_msets_smiles.txt", "w") as f:
        with Pool(32) as p:
            for file in p.imap_unordered(mset_to_smiles, glob.glob("/mnt/sdb1/adriscoll/ani1x-data/ani1x-msets/*.mset"), 100):
                if file is not None:
                    f.write(file+"\n")
