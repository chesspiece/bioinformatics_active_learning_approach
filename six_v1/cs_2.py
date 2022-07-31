import numpy as np
from numba.typed import List

from six_v1.comb_alg import cycle_to_chromosome

if __name__ == "__main__":
    with open("input_6_cs2.txt") as f:
        chromosome_graph = [
            int(x.strip("+")) for x in f.readline().strip("\n").strip("() ").split()
        ]

    chromosome = cycle_to_chromosome(List(chromosome_graph))

    with open("output.txt", "w") as f:
        #f.write("(" + " ".join(map(str, nodes)) + ")")
        f.write(" ".join(map(lambda x: f"{x:+d}", chromosome)))
