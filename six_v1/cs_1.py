import numpy as np
from numba.typed import List

from six_v1.comb_alg import chromosome_to_cycle

if __name__ == "__main__":
    with open("input_6_cs1.txt") as f:
        chromosome = [
            int(x.strip("+")) for x in f.readline().strip("\n").strip("() ").split()
        ]

    nodes = chromosome_to_cycle(List(chromosome))

    with open("output.txt", "w") as f:
        #f.write("(" + " ".join(map(str, nodes)) + ")")
        f.write(" ".join(map(str, nodes)))
