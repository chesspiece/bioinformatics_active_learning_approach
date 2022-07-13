import numpy as np
from numba.typed import List

from two.motif_lib import greedy_1

if __name__ == "__main__":
    with open("input.txt", "r") as f:
        dna_str = List(f.readline().strip())
        k = int(f.readline().strip())
        profile_mat = np.zeros((4, k), dtype=np.float64)

        for idx in range(4):
            profile_mat[idx, :] = np.array(
                [float(x) for x in f.readline().strip().split()]
            )
    res = greedy_1(profile_mat, dna_str, k)

    with open("output.txt", "w") as f:
        f.write("".join(res))
