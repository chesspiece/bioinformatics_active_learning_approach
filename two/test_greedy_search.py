import numpy as np
from numba.typed import List

from two.motif_lib import greedy_motif_search

if __name__ == "__main__":
    with open("greedy_motif_search_input.txt", "r") as f:
        k, t = (int(x) for x in f.readline().strip().split())
        dna_str = List(f.readline().strip().split())
    res = greedy_motif_search(dna_str, k, t)

    with open("output.txt", "w") as f:
        f.write(" ".join(res))
