import numpy as np

from six_v1.comb_alg import colored_edges, different_colored_edges

if __name__ == "__main__":
    with open("input_6_cs3.txt") as f:
        genome = f.readline().strip("\n").strip("()").split(")(")
        genome = [[int(synteny.strip("+")) for synteny in chromosome.split()] for chromosome in genome]
        

    edges = different_colored_edges(genome)
    edges = map(lambda x: (x[0], x[1]), edges)

    with open("output.txt", "w") as f:
        #f.write("(" + " ".join(map(str, nodes)) + ")")
        tmp = " ".join(map(lambda x: f"{x},", edges)).strip(",")
        f.write(tmp)
