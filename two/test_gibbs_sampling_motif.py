import numpy as np
from numba.typed import List

from two.motif_lib import gibbs_sampling_motif_finder

from tqdm import tqdm

if __name__ == "__main__":
    with open("gibbs_motif_find_input.txt", "r") as f:
        k, t, N = (int(x) for x in f.readline().strip().split())
        dna_str = List(f.readline().strip().split())
    best_res, best_score = gibbs_sampling_motif_finder(dna_str, k, N)
    for _ in tqdm(range(20)):
        res, score = gibbs_sampling_motif_finder(dna_str, k, N)
        if score < best_score:
            best_res = res
            best_score = score

    with open("output.txt", "w") as f:
        f.write(" ".join(best_res))
