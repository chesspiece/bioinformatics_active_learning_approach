import numpy as np
from six_v1.comb_alg import greedy_sorting

if __name__ == "__main__":
    with open("input_6_1.txt") as f:
        reversal = np.array(list(map(int, f.readline().strip().split())), dtype=np.int64)

    path = []
    dist, reversal_seq = greedy_sorting(reversal)

    with open("output.txt", "w") as f:
        f.write("\n".join([" ".join([f"{x:+d}" for x in reversal]) for reversal in reversal_seq]))
