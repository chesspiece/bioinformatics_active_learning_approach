from collections import defaultdict
import sys

from seven.tree_lib import small_parsimony

if __name__ == "__main__":
    with open("./input/input_7_6.txt") as f:
        n_leaves = int(f.readline().strip())
        tree: dict[str, list[str]] = defaultdict(list)
        leafs: dict[str, str] = defaultdict(str)
        for idx in range(n_leaves):
            first_node, second_node = f.readline().strip().split("->")
            tree[first_node].append(second_node)
            leafs[second_node] = second_node
        for idx, line in enumerate(f):
            first_node, second_node = line.strip().split("->")
            tree[first_node].append(second_node)
    try:
        root = first_node
    except NameError:
        print("Input file was empty")
        sys.exit()
    res = small_parsimony(tree, leafs, root)

