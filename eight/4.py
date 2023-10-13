import numpy as np

from clustering import hierarchical_clustering

if __name__ == "__main__":
    with open("input/input_8_4.txt", "r") as f:
        n = int(f.readline().strip())
        distance_matrix = np.zeros((n, n), dtype=np.float64)
        for idx, line in enumerate(f):
            line_data = list(map(float, line.strip().split()))
            distance_matrix[idx, :] = np.array(line_data)

    with open("output.txt", "w") as f:
        hierarchical_clustering(distance_matrix, f)
