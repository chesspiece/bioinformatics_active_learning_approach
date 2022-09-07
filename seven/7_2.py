import numpy as np

from seven.tree_lib import limb_length

if __name__ == "__main__":
    with open("./input/input_7_2.txt") as f:
        n_leaves = int(f.readline().strip())
        start_leafe = int(f.readline().strip())
        distance_matrix = np.zeros((n_leaves, n_leaves), dtype=np.int64)
        for idx, line in enumerate(f):
            distance_matrix[idx, :] = [int(x) for x in line.strip().split()]

    limb_len = limb_length(distance_matrix, start_leafe)
    with open("output.txt", "w") as f:
        f.write(str(limb_len))
