import glob
import json
import numpy as np
from rdkit import Chem, DataStructs

opt_smiles = []
try:
    with open("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/all_opt_smiles.txt", "r") as f:
        opt_smiles = json.loads(f.read())
except FileNotFoundError:
    for i in range(1, 10):
        for mol in glob.glob("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/opt/"+str(i)+"/*.mset"):
            mol_data = json.load(open(mol, "r"))
            try:
                opt_smiles.append(mol_data['identifiers']['smiles'])
            except KeyError:
                print(mol_data['identifiers'].keys(), mol)
                pass
    for x in ['b', 'br', 'i', 'p', 'p_new', 'se', 'si']:
        for mol in glob.glob("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/opt_cs40/"+x+"/*.mset"):
            mol_data = json.load(open(mol, "r"))
            try:
                opt_smiles.append(mol_data['identifiers']['smiles'])
            except KeyError:
                print(mol_data['identifiers'].keys(), mol)
                pass
    with open("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/all_opt_smiles.txt", "w") as file:
        json.dump(opt_smiles, file)

try:
    with open("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/all_opt_scores.npy", "rb") as f:
        opt_scores = np.load(f)
except FileNotFoundError:
    opt_mols = [Chem.MolFromSmiles(smile) for smile in opt_smiles]
    opt_fps = [Chem.RDKFingerprint(mol) for mol in opt_mols]
    opt_scores = np.zeros([len(opt_mols), len(opt_mols)])
    for i in range(len(opt_mols)):
        for j in range(i+1, len(opt_mols)):
            score = DataStructs.FingerprintSimilarity(opt_fps[i], opt_fps[j], metric=DataStructs.TanimotoSimilarity)
            opt_scores[i, j] = score
            opt_scores[j, i] = score
    with open("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/all_opt_scores.npy", "wb") as file:
        np.save(file, opt_scores)

meta_smiles = []
try:
    with open("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/all_meta_smiles.txt", "r") as f:
        meta_smiles = json.loads(f.read())
except FileNotFoundError:
    for mol in glob.glob("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/meta/*.mset"):
        mol_data = json.load(open(mol, "r"))
        try:
            meta_smiles.append(mol_data['identifiers']['smiles'])
        except KeyError:
            print(mol_data['identifiers'].keys(), mol)
            pass
    for mol in glob.glob("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/aimd/*.mset"):
        mol_data = json.load(open(mol, "r"))
        try:
            meta_smiles.append(mol_data['identifiers']['smiles'])
        except KeyError:
            print(mol_data['identifiers'].keys(), mol)
            pass
    with open("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/all_meta_smiles.txt", "w") as file:
        json.dump(meta_smiles, file)


try:
    with open("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/all_meta_scores.npy", "rb") as f:
        meta_scores = np.load(f)
except FileNotFoundError:
    meta_mols = [Chem.MolFromSmiles(smile) for smile in meta_smiles]
    meta_fps = [Chem.RDKFingerprint(mol) for mol in meta_mols]
    meta_scores = np.zeros([len(meta_mols), len(meta_mols)])
    for i in range(len(meta_mols)):
        for j in range(i+1, len(meta_mols)):
            score = DataStructs.FingerprintSimilarity(meta_fps[i], meta_fps[j], metric=DataStructs.TanimotoSimilarity)
            meta_scores[i, j] = score
            meta_scores[j, i] = score
    with open("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/all_meta_scores.npy", "wb") as file:
        np.save(file, meta_scores)

ani_smiles = []
try:
    with open("/mnt/sdb1/adriscoll/ani1x-data/all_ani_smiles.txt", "r") as f:
        ani_smiles = json.loads(f.read())
except FileNotFoundError:
    for mol in glob.glob("/mnt/sdb1/adriscoll/ani1x-data/ani1x-msets/*.mset"):
        mol_data = json.load(open(mol, "r"))
        try:
            ani_smiles.append(mol_data['identifiers']['smiles'])
        except KeyError:
            print(mol_data['identifiers'].keys(), mol)
            pass
    with open("/mnt/sdb1/adriscoll/ani1x-data/all_ani_smiles.txt", "w") as file:
        json.dump(ani_smiles, file)

try:
    with open("/mnt/sdb1/adriscoll/ani1x-data/all_ani_scores.npy", "rb") as f:
        ani_scores = np.load(f)
except FileNotFoundError:
    ani_mols = [Chem.MolFromSmiles(smile) for smile in ani_smiles]
    ani_fps = [Chem.RDKFingerprint(mol) for mol in ani_mols]
    ani_scores = np.zeros([len(ani_mols), len(ani_mols)])
    for i in range(len(ani_mols)):
        for j in range(i+1, len(ani_mols)):
            score = DataStructs.FingerprintSimilarity(ani_fps[i], ani_fps[j], metric=DataStructs.TanimotoSimilarity)
            ani_scores[i, j] = score
            ani_scores[j, i] = score
    with open("/mnt/sdb1/adriscoll/ani1x-data/all_ani_scores.npy", "wb") as file:
        np.save(file, ani_scores)

gdb9_smiles = []
try:
    with open("/mnt/sdb1/adriscoll/gdb9-data/all_gdb9_smiles.txt", "r") as f:
        gdb9_smiles = json.loads(f.read())
except FileNotFoundError:
    for mol in glob.glob("/mnt/sdb1/adriscoll/gdb9-data/dsgdb9nsd/*.xyz"):
        mol_data = (open(mol, "r")).readlines()
        gdb9_smiles.append(mol_data[-2].split()[0])
    with open("/mnt/sdb1/adriscoll/gdb9-data/all_gdb9_smiles.txt", "w") as file:
        json.dump(gdb9_smiles, file)

try:
    with open("/mnt/sdb1/adriscoll/gdb9-data/all_gdb9_scores.npy", "rb") as f:
        gdb9_scores = np.load(f)
except FileNotFoundError:
    gdb9_mols = [Chem.MolFromSmiles(smile) for smile in gdb9_smiles]
    gdb9_fps = [Chem.RDKFingerprint(mol) for mol in gdb9_mols]
    gdb9_scores = np.zeros([len(gdb9_mols), len(gdb9_mols)])
    for i in range(len(gdb9_mols)):
        for j in range(i+1, len(gdb9_mols)):
            score = DataStructs.FingerprintSimilarity(gdb9_fps[i], gdb9_fps[j], metric=DataStructs.TanimotoSimilarity)
            gdb9_scores[i, j] = score
            gdb9_scores[j, i] = score
    with open("/mnt/sdb1/adriscoll/gdb9-data/all_gdb9_scores.npy", "wb") as file:
        np.save(file, gdb9_scores)

