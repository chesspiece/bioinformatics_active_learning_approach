import numpy as np
from six_v1.comb_alg import breakpoints_number

if __name__ == "__main__":
    with open("input_6_2.txt") as f:
        reversal = np.array(list(map(int, f.readline().strip().split())), dtype=np.int64)

    path = []
    breakpoints_quant = breakpoints_number(reversal)

    with open("output.txt", "w") as f:
        f.write(str(breakpoints_quant))
