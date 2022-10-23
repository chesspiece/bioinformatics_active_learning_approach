import sys
from collections import defaultdict

from one.skew import hamming_str
from seven.tree_lib import small_parsimony, undirected2directed, insert_root

if __name__ == "__main__":
    with open("./input/input_7_7.txt", "r") as f:
        n_leaves = int(f.readline().strip())
        tree: dict[str, list[str]] = defaultdict(list)
        leafs: dict[str, str] = defaultdict(str)
        for idx in range(2 * n_leaves):
            # f.readline().strip().split("->")
            first_node, second_node = f.readline().strip().split("->")
            tree[first_node].append(second_node)
            if (
                "A" in second_node
                or "C" in second_node
                or "G" in second_node
                or "T" in second_node
            ):
                leafs[second_node] = second_node
        for idx, line in enumerate(f):
            # if idx % 2 == 0:
            #    continue
            first_node, second_node = line.strip().split("->")
            tree[first_node].append(second_node)
    try:
        root = insert_root(tree)
    except NameError:
        print("Input file was empty")
        sys.exit()
    undirected2directed(tree, root)
    final_score, true_tree = small_parsimony(tree, leafs, root, additional_root=True)
    root_son, root_daughter = true_tree[root]
    true_tree[root_son].append(root_daughter)
    true_tree[root_daughter].append(root_son)
    with open("output.txt", "w") as f:
        f.write(f"{final_score}\n")
        for node, node_neighbours in true_tree.items():
            for neighbour in node_neighbours:
                if node == root or neighbour == root:
                    continue
                f.write(f"{node}->{neighbour}:{hamming_str(node, neighbour)}\n")
