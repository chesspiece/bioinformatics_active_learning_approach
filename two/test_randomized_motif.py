import numpy as np
from numba.typed import List

from two.motif_lib import randomized_motif_search

from tqdm import tqdm

if __name__ == "__main__":
    with open("randomized_motif_search_input.txt", "r") as f:
        k, t = (int(x) for x in f.readline().strip().split())
        dna_str = List(f.readline().strip().split())
    best_res, best_score = randomized_motif_search(dna_str, k)
    for _ in tqdm(range(729)):
        res, score = randomized_motif_search(dna_str, k)
        if score < best_score:
            best_res = res
            best_score = score

    with open("output.txt", "w") as f:
        f.write(" ".join(best_res))
