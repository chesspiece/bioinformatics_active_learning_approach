import numpy as np

from five.dp_lib import manhattan_tourist

if __name__ == "__main__":
    with open("input_5_2.txt") as f:
        n, m = (int(x) for x in f.readline().strip().split())
        down = np.zeros((n, m + 1), dtype=np.int64)
        right = np.zeros((n + 1, m), dtype=np.int64)
        for i in range(n):
            down[i, :] = [int(x) for x in f.readline().strip().split()]
        f.readline()
        for i in range(n + 1):
            right[i, :] = [int(x) for x in f.readline().strip().split()]
    n += 1
    m += 1

    answ = manhattan_tourist(n, m, down, right)

    with open("output.txt", "w") as f:
        f.write("".join(str(answ)))
