import numpy as np

from six_v1.comb_alg import two_breaks_genome

if __name__ == "__main__":
    with open("input_6_cs6.txt") as f:
        genome = f.readline().strip("\n").strip("()").split(")(")
        genome = [[int(synteny.strip("+")) for synteny in chromosome.split()] for chromosome in genome]
        i, j, k, l = [int(x) for x in f.readline().strip().split(", ")]
        

    genome = two_breaks_genome(genome, i, j, k, l)
    genome = [" ".join([f"{gen:+d}" for gen in chromosome]) for chromosome in genome]

    with open("output.txt", "w") as f:
        # f.write("(" + " ".join(map(str, nodes)) + ")")
        tmp = "(" + ")(".join(genome) + ")"
        f.write(tmp)
