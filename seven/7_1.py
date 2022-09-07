from collections import defaultdict

import numpy as np

from seven.tree_lib import dijkstra

if __name__ == "__main__":
    with open("./input/input_7_1.txt") as f:
        n_leaves = int(f.readline().strip())
        adj_list: dict[int, dict[int, int]] = defaultdict(dict)
        node_quant = n_leaves - 1
        for line in f:
            start_node, end_node = line.strip().split("->")
            start_node = int(start_node)
            end_node, weight = (int(x) for x in end_node.split(":"))
            if start_node > node_quant:
                node_quant = start_node
            if (end_node, weight) in adj_list[start_node].items():
                continue
            else:
                adj_list[start_node][end_node] = weight
                adj_list[end_node][start_node] = weight

    with open("output.txt", "w") as f:
        for i in range(n_leaves):
            dists, _ = dijkstra(adj_list, i, node_quant + 1)
            dists = dists[0:n_leaves]
            f.write(" ".join(map(str, dists)))
            f.write("\n")
