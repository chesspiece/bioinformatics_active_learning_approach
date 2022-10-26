import sys
from collections import defaultdict

from one.skew import hamming_str
from seven.tree_lib import nearest_neighbour_tree

if __name__ == "__main__":
    with open("./input/input_7_8.txt", "r") as f:
        edge_1, edge_2 = f.readline().strip().split()
        tree: dict[str, list[str]] = defaultdict(list)
        for idx, line in enumerate(f):
            first_node, second_node = line.strip().split("->")
            tree[first_node].append(second_node)
    tree_1, tree_2 = nearest_neighbour_tree(tree, edge_1, edge_2)
    with open("output.txt", "w") as f:
        for node, node_neighbours in tree_1.items():
            for neighbour in node_neighbours:
                f.write(f"{node}->{neighbour}\n")
        f.write("\n")
        for node, node_neighbours in tree_2.items():
            for neighbour in node_neighbours:
                f.write(f"{node}->{neighbour}\n")
