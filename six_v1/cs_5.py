from six_v1.comb_alg import two_breaks_genome_graph

if __name__ == "__main__":
    with open("input_6_cs5.txt") as f:
        genome_edges = f.readline().strip("\n").strip("()").split("), (")
        genome_edges = [x.split(", ") for x in genome_edges]
        genome_edges = [(int(x), int(y)) for x, y in genome_edges]
        i, j, k, l = [int(x) for x in f.readline().strip().split(", ")]

    genome_edges = two_breaks_genome_graph(genome_edges, i, j, k, l)

    with open("output.txt", "w") as f:
        tmp = " ".join(map(lambda x: f"{x},", genome_edges)).strip(",")
        f.write(tmp)
