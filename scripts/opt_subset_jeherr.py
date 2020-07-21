import json
import time
import numpy as np
from itertools import combinations
from rdkit import Chem, DataStructs
from multiprocessing import Pool


def build_similarity_matrix():
    with open("/home/jeherr/tensorchem/tmp/all_opt_smiles.txt", "r") as f:
        lines = f.readlines()
        opt_smiles = [line.strip("\n") for line in lines]

    opt_mols = [Chem.MolFromSmiles(smile) for smile in opt_smiles]
    rarest_elements = [5, 14, 15, 34, 35, 53]
    keep_idx = [i for i, mol in enumerate(opt_mols) if any(item in rarest_elements for item in mol.GetAtoms())]
    #opt_less = [mol for mol in opt_mols if mol.GetNumHeavyAtoms() < 25]
    opt_fps = [Chem.RDKFingerprint(mol) for mol in opt_mols]

    opt_scores = np.ones((len(opt_mols), len(opt_mols)))
    for i, fp1 in enumerate(opt_fps):
        for j, fp2 in enumerate(opt_fps[i+1:]):
            score = DataStructs.FingerprintSimilarity(fp1, fp2, metric=DataStructs.TanimotoSimilarity)
            opt_scores[i,i+j+1] = opt_scores[i+j+1,i] = score

    np.save("/home/jeherr/tensorchem/tmp/all_opt_scores.npy", opt_scores)
    opt_scores = np.square(opt_scores)
    np.save("/home/jeherr/tensorchem/tmp/opt_scores_squared.npy", opt_scores)
    return keep_idx


keep_idx = build_similarity_matrix()


def optimize_subset(k, keep_idx):
    opt_scores = np.load("/home/jeherr/tensorchem/tmp/opt_scores_squared.npy")
    opt_scores_sum = np.sum(opt_scores, axis=1)
    min_idx = np.argsort(opt_scores_sum).tolist()
    min_idx_less_keep = [idx for idx in min_idx if idx not in keep_idx]
    keep_len = len(keep_idx)
    additional_len = k - keep_len
    min_k = min_idx_less_keep[:additional_len]
    sorted_idx = min_idx_less_keep[additional_len:]

    keep_combos = combinations(keep_idx, 2)
    keep_score = sum([opt_scores[combo[0], combo[1]] for combo in keep_combos])

    best_combos = combinations(min_k, 2)
    init_score = 0
    for combo in best_combos:
        init_score += opt_scores[combo[0], combo[1]]

    best_score = init_score

    while True:
        none_changed = True
        for i, orig_idx in enumerate(min_k[keep_len:]):
            unchanged_idx = [idx for idx in min_k if idx != orig_idx]
            old_combos = [[orig_idx, idx] for idx in unchanged_idx]
            old_score = sum([opt_scores[combo[0], combo[1]] for combo in old_combos])
            for j, new_idx in enumerate(sorted_idx):
                new_score = sum([opt_scores[new_idx, const_idx] for const_idx in unchanged_idx])
                new_best = best_score - old_score + new_score
                if new_best < best_score:
                    print(new_best / k, best_score / k)
                    none_changed = False
                    best_score = new_best
                    old_score = new_score
                    min_k[keep_len+i] = new_idx
                    sorted_idx[j] = orig_idx
        if none_changed:
            break

    opt_scores = np.load("/home/jeherr/tensorchem/tmp/all_opt_scores.npy")
    subopt_scores = np.ones((k, k))
    for i, idx1 in enumerate(min_k):
        for j, idx2 in enumerate(min_k[i+1:]):
            subopt_scores[i,j] = subopt_scores[j,i] = opt_scores[idx1,idx2]

    np.save(f"/home/jeherr/tensorchem/tmp/subopt_scores_{k}_{init_score}_{best_score}.npy", subopt_scores)
    with open(f"/home/jeherr/tensorchem/tmp/subopt_scores_{k}_{init_score}_{best_score}.txt", "w") as f:
        for idx in min_k:
            f.write(idx)


optimize_subset(6000, keep_idx)
