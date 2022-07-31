import numpy as np

from six_v1.comb_alg import shortest_rearrangmen

if __name__ == "__main__":
    with open("input_6_4.txt") as f:
        genome_1 = f.readline().strip("\n").strip("()").split(")(")
        genome_1 = [[int(synteny.strip("+")) for synteny in chromosome.split()] for chromosome in genome_1]
        genome_2 = f.readline().strip("\n").strip("()").split(")(")
        genome_2 = [[int(synteny.strip("+")) for synteny in chromosome.split()] for chromosome in genome_2]
        
    res = shortest_rearrangmen(genome_1, genome_2)
    with open("output.txt", "w") as f:
        f.write("\n".join(res))
