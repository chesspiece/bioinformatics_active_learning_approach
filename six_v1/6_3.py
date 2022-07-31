import numpy as np

from six_v1.comb_alg import colored_edges, different_colored_edges, create_adjacency_list, dfs_find_cicle_black

if __name__ == "__main__":
    with open("input_6_3.txt") as f:
        genome_1 = f.readline().strip("\n").strip("()").split(")(")
        genome_1 = [[int(synteny.strip("+")) for synteny in chromosome.split()] for chromosome in genome_1]
        genome_2 = f.readline().strip("\n").strip("()").split(")(")
        genome_2 = [[int(synteny.strip("+")) for synteny in chromosome.split()] for chromosome in genome_2]
        
    edges_1 = different_colored_edges(genome_1, color=0)
    edges_2 = different_colored_edges(genome_2, color=1)
    edges = edges_1 + edges_2
    adj_list = create_adjacency_list(edges)
    cycles = dfs_find_cicle_black(adj_list)

    with open("output.txt", "w") as f:
        f.write(str(len(edges_1) - len(cycles)))
