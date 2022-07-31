import numpy as np

from numba.typed import List
from six_v1.comb_alg import edges_to_genome

if __name__ == "__main__":
    with open("input_6_cs4.txt") as f:
        genome_edges = f.readline().strip("\n").strip("()").split("), (")
        genome_edges = [x.split(", ") for x in genome_edges]
        genome_edges = [(int(x), int(y)) for x, y in genome_edges]

    genome = edges_to_genome(genome_edges)
    genome = [" ".join([f"{gen:+d}" for gen in chromosome]) for chromosome in genome]

    with open("output.txt", "w") as f:
        # f.write("(" + " ".join(map(str, nodes)) + ")")
        tmp = "(" + ")(".join(genome) + ")"
        f.write(tmp)
