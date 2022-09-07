import numpy as np

from seven.tree_lib import additive_philogeny

if __name__ == "__main__":
    with open("./input/input_7_3.txt") as f:
        n_leaves = int(f.readline().strip())
        distance_matrix = np.zeros((n_leaves, n_leaves), dtype=np.int64)
        for idx, line in enumerate(f):
            distance_matrix[idx, :] = [int(x) for x in line.strip().split()]

    trr, num_nodes = additive_philogeny(distance_matrix, n_leaves)
    with open("output.txt", "w") as f:
        for i in range(num_nodes):
            for key, val in trr[i].items():
                f.write(f"{i}->{key}:{int(val)}\n")
