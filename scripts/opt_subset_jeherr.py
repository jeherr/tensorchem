import json
import time
import numpy as np
from itertools import combinations
from rdkit import Chem, DataStructs
from multiprocessing import Pool


def build_similarity_matrix():
    with open("/home/jeherr/tensorchem/tmp/all_opt_smiles.txt", "r") as f:
        opt_smiles = json.loads(f.read())
        #sub_opt_smiles = opt_smiles[:10000]

    opt_mols = [Chem.MolFromSmiles(smile) for smile in opt_smiles]
    opt_less = [mol for mol in opt_mols if mol.GetNumHeavyAtoms() < 25]
    opt_fps = [Chem.RDKFingerprint(mol) for mol in opt_less]

    opt_scores = np.ones((len(opt_less), len(opt_less)))
    for i, fp1 in enumerate(opt_fps):
        for j, fp2 in enumerate(opt_fps[i+1:]):
            score = DataStructs.FingerprintSimilarity(fp1, fp2, metric=DataStructs.TanimotoSimilarity)
            opt_scores[i,i+j+1] = opt_scores[i+j+1,i] = score

    np.save("/home/jeherr/tensorchem/tmp/all_opt_scores.npy", opt_scores)
    opt_scores = np.square(opt_scores)
    np.save("/home/jeherr/tensorchem/tmp/opt_scores_squared.npy", opt_scores)


#build_similarity_matrix()


def optimize_subset(k):
    opt_scores = np.load("/home/jeherr/tensorchem/tmp/opt_scores_squared.npy")
    opt_scores_sum = np.sum(opt_scores, axis=1)
    min_idx = np.argsort(opt_scores_sum)
    min_k = min_idx[:k].tolist()
    sorted_idx = min_idx[k:].tolist()

    best_combos = combinations(min_k, 2)
    init_score = 0
    for combo in best_combos:
        init_score += opt_scores[combo[0], combo[1]]

    best_score = init_score

    while True:
        none_changed = True
        for i, orig_idx in enumerate(min_k):
            unchanged_idx = [idx for idx in min_k if idx != orig_idx]
            old_combos = [[orig_idx, idx] for idx in unchanged_idx]
            # unchanged_combos = combinations(unchanged_idx, 2)
            old_score = sum([opt_scores[combo[0], combo[1]] for combo in old_combos])
            for j, new_idx in enumerate(sorted_idx):
                new_score = sum([opt_scores[new_idx, const_idx] for const_idx in unchanged_idx])
                new_best = best_score - old_score + new_score
                if new_best < best_score:
                    #print(new_best / k, best_score / k)
                    none_changed = False
                    best_score = new_best
                    old_score = new_score
                    min_k[i] = new_idx
                    sorted_idx[j] = orig_idx
        if none_changed:
            break

    opt_scores = np.load("/home/jeherr/tensorchem/tmp/all_opt_scores.npy")
    subopt_scores = np.ones((k, k))
    for i, idx1 in enumerate(min_k):
        for j, idx2 in enumerate(min_k[i+1:]):
            subopt_scores[i,j] = subopt_scores[j,i] = opt_scores[idx1,idx2]

    np.save(f"/home/jeherr/tensorchem/tmp/subopt_scores_{k}_{init_score}_{best_score}.npy", subopt_scores)


with Pool(20) as p:
    p.map(optimize_subset, [250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000, 3250, 3500, 3750, 4000, 4250, 4500, 4750, 5000])
