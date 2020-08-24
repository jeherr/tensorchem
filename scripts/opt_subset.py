import json
import time
import numpy as np
from itertools import combinations
from rdkit import Chem, DataStructs

with open("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/all_opt_smiles.txt", "r") as f:
    opt_smiles = json.loads(f.read())

opt_mols = [Chem.MolFromSmiles(smile) for smile in opt_smiles]
opt_fps = [Chem.RDKFingerprint(mol) for mol in opt_mols]

with open("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/all_opt_scores.npy", "rb") as f:
    opt_scores = np.load(f)

opt_scores = np.square(opt_scores)

with open("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/opt_init_indexes.txt", "r") as f:
    init_indexes = json.loads(f.read())

k = 3000
init_sum_score = 0
for i in range(len(init_indexes)):
    row_score = 0
    for j in range(i+1, len(init_indexes)):
        if init_indexes[i] == init_indexes[j]:
            row_score += 1
        else:
            score = DataStructs.FingerprintSimilarity(opt_fps[init_indexes[i]], opt_fps[init_indexes[j]], metric=DataStructs.TanimotoSimilarity)
            row_score += 2 * score * score
    init_sum_score += row_score

print("initial min score: ", init_sum_score)
min_sum_score = init_sum_score

all_idx = init_indexes.copy()
print("done...starting list...")

for x in range(55804):
    if x not in init_indexes:
        sum_scores = []
        for idx in range(len(all_idx)):
            temp = all_idx[idx]
            all_idx[idx] = x
            sum_score = 0
            for i in range(len(all_idx)):
                for j in range(i+1, len(all_idx)):
                    if all_idx[i] == all_idx[j]:
                        sum_score += 1
                    else:
                        sum_score += 2 * opt_scores[all_idx[i]][all_idx[j]]
            all_idx[i] = temp
        print(min(sum_scores), x, sum_scores.index(min_sum_score))
        if min(sum_scores) < min_sum_score:
            min_sum_score = min(sum_scores)
            all_idx[sum_scores.index(min_sum_score)] = x
            print("new min score! using ", x, " at index ", sum_scores.index(min_sum_score))

min_indexes = all_idx.copy()
opt_sub_scores = np.zeros([k, k])
for i in range(len(min_indexes)):
    for j in range(i+1, len(min_indexes)):
        score = DataStructs.FingerprintSimilarity(opt_fps[min_indexes[i]], opt_fps[min_indexes[j]], metric=DataStructs.TanimotoSimilarity)
        opt_sub_scores[i, j] = score
        opt_sub_scores[j, i] = score

print("final min score: ", np.sum(np.square(opt_sub_scores)))

with open("/mnt/sdb1/adriscoll/chemspider_data/expanded_msets/opt_subset_"+str(k)+"_scores.npy", "wb") as f:
    np.save(f, opt_sub_scores)

