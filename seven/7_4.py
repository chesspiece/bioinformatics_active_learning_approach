import numpy as np

from seven.tree_lib import ugma_philogeny

if __name__ == "__main__":
    with open("./input/input_7_4.txt") as f:
        n_leaves = int(f.readline().strip())
        distance_matrix = np.zeros((n_leaves, n_leaves), dtype=np.float64)
        for idx, line in enumerate(f):
            distance_matrix[idx, :] = [int(x) for x in line.strip().split()]

    trr = ugma_philogeny(distance_matrix)
    with open("output.txt", "w") as f:
        for k in sorted(trr.keys()):
            for val, age_dur in trr[k].items():
                f.write(f"{k}->{val}:{age_dur:.3f}\n")
