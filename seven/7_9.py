import sys
from collections import defaultdict

from one.skew import hamming_str
from seven.tree_lib import (delete_root, insert_root,
                            nearest_neighbor_parsimony, small_parsimony,
                            undirected2directed)

if __name__ == "__main__":
    with open("./input/input_7_9.txt", "r") as f:
        n_leaves = int(f.readline().strip())
        tree: dict[str, list[str]] = defaultdict(list)
        leafs: dict[str, str] = defaultdict(str)
        for idx, line in enumerate(f):
            # f.readline().strip().split("->")
            first_node, second_node = line.strip().split("->")
            tree[first_node].append(second_node)
            if (
                "A" in second_node
                or "C" in second_node
                or "G" in second_node
                or "T" in second_node
            ):
                leafs[second_node] = second_node
    with open("output.txt", "w") as f:
        nearest_neighbor_parsimony(tree, leafs, f)
