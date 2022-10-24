from collections import defaultdict
from copy import deepcopy
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import scipy as scp

from mapped_query import MappedQueue


def dijkstra(adj_list: dict[int, dict[int, int]], start_edge: int, v_quant: int):
    dist = [np.inf for _ in range(v_quant)]
    prev = [-1 for _ in range(v_quant)]
    dist[start_edge] = 0
    dist_heap = [(dist[i], i) for i in adj_list.keys()]
    dist_heap = MappedQueue(dist_heap)
    while dist_heap:
        _, vertex = dist_heap.pop()
        for u, weight in adj_list[vertex].items():
            alt = dist[vertex] + weight
            if alt < dist[u]:
                updated_elem = (dist[u], u)
                dist[u] = alt
                prev[u] = vertex
                dist_heap.update(updated_elem, (alt, u))
    return dist, prev


def get_path_sequence(prevs: list[int], i, j) -> list[int]:
    path = [j]
    prev_node = prevs[j]
    while prev_node != i and prev_node != -1:
        path.append(prev_node)
        prev_node = prevs[prev_node]
    if prev_node == -1:
        raise ValueError("Error")
    path.append(i)
    return path[::-1]


def limb_length(distance_matrix: npt.NDArray, start_leafe: int):
    n_leafes, _ = distance_matrix.shape
    min_dist = np.inf
    j = 0 if start_leafe != 0 else 1
    for i in range(n_leafes):
        # for j in range(n_leafes):
        if i == start_leafe or j == start_leafe or i == j:
            continue
        test_dist = (
            distance_matrix[i, start_leafe]
            + distance_matrix[start_leafe, j]
            - distance_matrix[i, j]
        )
        if test_dist < min_dist:
            min_dist = test_dist
    return min_dist // 2


def find_path_nodes(matrix, n, limb_len):
    """
    Finds two nodes (i, j) such that they are located in two different subtrees that is joined by the parent node of n.

    :param matrix: A distance matrix
    :param n: Node dividing the two subtrees
    :param limb_len: Limb length of node n
    :return: Two nodes in the graph
    """
    limb_len *= 2
    for i in range(n - 1):
        for j in range(n - 1):
            if matrix[i][j] == matrix[i][n - 1] + matrix[n - 1][j] - limb_len:
                return i, j
    raise ValueError("Wrong matrix!")


def additive_philogeny(
    distane_matrix: npt.NDArray[np.int64],
    num_nodes: int,
) -> tuple[dict[int, dict[int, int]], int]:
    n, _ = distane_matrix.shape
    if n == 2:
        ret = {0: {1: distane_matrix[0, 1]}, 1: {0: distane_matrix[1, 0]}}
        return ret, num_nodes

    tree, num_nodes = additive_philogeny(distane_matrix[:-1, :-1], num_nodes)

    limb_lngth = int(limb_length(distane_matrix, n - 1))
    i, j = find_path_nodes(distane_matrix, n, limb_lngth)

    add_vertex_dist = distane_matrix[i, n - 1] - limb_lngth

    path_length, all_prevs = dijkstra(tree, i, num_nodes)
    path = get_path_sequence(all_prevs, i, j)
    curr_dist = 0
    for node_idx in range(len(path) - 1):
        curr_node = path[node_idx]
        next_node = path[node_idx + 1]
        curr_dist = path_length[next_node]
        if curr_dist > add_vertex_dist:
            tree[next_node].pop(curr_node)
            edge_len = tree[curr_node].pop(next_node)

            tree[n - 1] = {num_nodes: limb_lngth}
            tree[num_nodes] = {n - 1: limb_lngth}

            tree[next_node][num_nodes] = curr_dist - add_vertex_dist
            tree[curr_node][num_nodes] = edge_len - (curr_dist - add_vertex_dist)

            tree[num_nodes][next_node] = curr_dist - add_vertex_dist
            tree[num_nodes][curr_node] = edge_len - (curr_dist - add_vertex_dist)

            num_nodes += 1
            break
        elif curr_dist == add_vertex_dist:
            tree[next_node][n] = limb_lngth
            tree[n - 1] = {next_node: limb_lngth}
            break
    return tree, num_nodes


def ugma_philogeny(
    distane_matrix: npt.NDArray[np.floating],
) -> dict[int, dict[int, int]]:
    max_node, _ = distane_matrix.shape
    clusters = {i: [i] for i in range(max_node)}
    tree: dict[int, dict[int, int]] = {i: {} for i in range(max_node)}
    age = {i: 0 for i in range(max_node)}

    np.fill_diagonal(distane_matrix, np.inf)

    while len(clusters) != 1:
        closest_node_i, closest_node_j = (
            int(x)
            for x in np.unravel_index(np.argmin(distane_matrix), distane_matrix.shape)
        )
        new_age = distane_matrix[closest_node_i, closest_node_j] / 2

        root_i = list(clusters.keys())[closest_node_i]
        root_j = list(clusters.keys())[closest_node_j]

        cluster_i = clusters.pop(root_i)
        cluster_j = clusters.pop(root_j)
        new_cluster = cluster_i + cluster_j
        clusters[max_node] = new_cluster

        new_cluster_dists = (
            len(cluster_i)
            * np.delete(
                distane_matrix[:, closest_node_i], [closest_node_i, closest_node_j]
            )
            + len(cluster_j)
            * np.delete(
                distane_matrix[:, closest_node_j], [closest_node_j, closest_node_i]
            )
        ) / (len(cluster_i) + len(cluster_j))

        age[max_node] = new_age
        tree[max_node] = {root_i: new_age - age[root_i], root_j: new_age - age[root_j]}
        tree[root_i][max_node] = new_age - age[root_i]
        tree[root_j][max_node] = new_age - age[root_j]

        distane_matrix = np.delete(
            distane_matrix, [closest_node_i, closest_node_j], axis=0
        )
        distane_matrix = np.delete(
            distane_matrix, [closest_node_i, closest_node_j], axis=1
        )
        distane_matrix = np.vstack((distane_matrix, new_cluster_dists.reshape(1, -1)))
        distane_matrix = np.hstack(
            (distane_matrix, np.append(new_cluster_dists, np.inf).reshape(-1, 1))
        )

        max_node += 1

    return tree


def neighbours_joining_phylogeny(
    distance_matrix: npt.NDArray[np.floating],
) -> dict[int, dict[int, int]]:
    max_nodes, _ = distance_matrix.shape
    nodes_dict = {i: 0 for i in range(max_nodes)}
    max_nodes -= 1
    tree, _ = __neighbours_joining_phylogeny__(distance_matrix, max_nodes, nodes_dict)
    return tree


def __neighbours_joining_phylogeny__(
    distance_matrix: npt.NDArray[np.floating], max_node: int, nodes_dict: dict[int, int]
) -> tuple[dict[int, dict[int, int]], int]:
    num_nodes, _ = distance_matrix.shape

    if num_nodes == 2:

        node_i = list(nodes_dict.keys())[0]
        node_j = list(nodes_dict.keys())[1]

        _ = nodes_dict.pop(node_i)
        _ = nodes_dict.pop(node_j)

        ret = {
            node_i: {node_j: distance_matrix[0, 1]},
            node_j: {node_i: distance_matrix[1, 0]},
        }
        return ret, max_node

    total_distance = distance_matrix.sum(axis=1)
    distance_matrix_star = (
        (num_nodes - 2) * distance_matrix
        - total_distance
        - total_distance.reshape(-1, 1)
    )
    np.fill_diagonal(distance_matrix_star, np.inf)

    closest_node_i, closest_node_j = (
        int(x)
        for x in np.unravel_index(
            np.argmin(distance_matrix_star), distance_matrix_star.shape
        )
    )
    # node_i = nodes_dict.pop(closest_node_i)
    # node_j = nodes_dict.pop(closest_node_j)

    node_i = list(nodes_dict.keys())[closest_node_i]
    node_j = list(nodes_dict.keys())[closest_node_j]

    _ = nodes_dict.pop(node_i)
    _ = nodes_dict.pop(node_j)

    delta = (total_distance[closest_node_i] - total_distance[closest_node_j]) / (
        num_nodes - 2
    )
    limb_length_i = 0.5 * (distance_matrix[closest_node_i, closest_node_j] + delta)
    limb_length_j = 0.5 * (distance_matrix[closest_node_i, closest_node_j] - delta)

    new_cluster_dists = 0.5 * (
        distance_matrix[:, closest_node_i]
        + distance_matrix[:, closest_node_j]
        - distance_matrix[closest_node_i, closest_node_j]
    )
    new_cluster_dists = np.delete(new_cluster_dists, [closest_node_i, closest_node_j])

    distance_matrix = np.delete(
        distance_matrix, [closest_node_i, closest_node_j], axis=0
    )
    distance_matrix = np.delete(
        distance_matrix, [closest_node_i, closest_node_j], axis=1
    )

    distance_matrix = np.vstack((distance_matrix, new_cluster_dists.reshape(1, -1)))
    distance_matrix = np.hstack(
        (distance_matrix, np.append(new_cluster_dists, 0).reshape(-1, 1))
    )

    max_node += 1
    nodes_dict[max_node] = 0
    tree, _ = __neighbours_joining_phylogeny__(distance_matrix, max_node, nodes_dict)

    # tree[max_node] = {node_i: limb_length_i, node_j: limb_length_j}
    tree[max_node][node_i] = limb_length_i
    tree[max_node][node_j] = limb_length_j
    # tree[node_i][max_node] = limb_length_i
    # tree[node_j][max_node] = limb_length_j
    tree[node_i] = {max_node: limb_length_i}
    tree[node_j] = {max_node: limb_length_j}
    return tree, max_node


def is_leaf(tree: dict[str, list[str]], node: str) -> bool:
    if not tree[node]:
        return True
    return False


def __backtrack__(
    tree: dict[str, list[str]],
    backtrack,
    root: str,
    letter,
    ddct: dict[str, str],
    char_idx: int,
) -> dict[str, str]:
    ncl_dict = {"A": 0, "C": 1, "G": 2, "T": 3}
    ncl_dict_rev = {ncl_dict[k]: k for k in ncl_dict.keys()}
    if is_leaf(tree, root):
        ddct[root] = root[char_idx]
        return ddct
    daughter, son = tree[root]
    letter1, letter2 = backtrack[root][0]
    __backtrack__(tree, backtrack, daughter, letter1[letter], ddct, char_idx)
    __backtrack__(tree, backtrack, son, letter2[letter], ddct, char_idx)
    ddct[root] = ncl_dict_rev[letter]
    return ddct


def __small_parsimony__(
    tree: dict[str, list[str]],
    leaf_characters: dict[str, str],
    root: str,
    char_idx: int,
) -> tuple[int, dict[str, str]]:
    ncl_dict = {"A": 0, "C": 1, "G": 2, "T": 3}

    leaf_keys = list(leaf_characters.keys())
    keys = list(tree.keys())
    tags = {k: 0 for k in leaf_keys + keys}
    score = {k: np.array([1] * 4) for k in leaf_keys + keys}
    backtrack = {k: [] for k in leaf_keys + keys}

    a = np.ones((4, 4))
    np.fill_diagonal(a, 0)

    for node in leaf_keys:
        tags[node] = 1
        score[node][ncl_dict[leaf_characters[node][char_idx]]] = 0

    for _ in range(len(tree.keys())):
        for node in keys:
            try:
                daughter, son = tree[node]
            except Exception:
                continue
            if tags[node] == 0 and tags[daughter] == 1 and tags[son] == 1:
                score[node] = np.min(score[son] + a, axis=1) + np.min(
                    score[daughter] + a, axis=1
                )
                backtrack[node].append(
                    [
                        np.argmin(score[daughter] + a, axis=1),
                        np.argmin(score[son] + a, axis=1),
                    ]
                )
                tags[node] = 1
    t_dct = defaultdict(str)
    res = __backtrack__(tree, backtrack, root, np.argmin(score[root]), t_dct, char_idx)
    return int(np.min(score[root])), res


def small_parsimony(
    tree: dict[str, list[str]],
    leaf_characters: dict[str, str],
    root: str,
    additional_root=False,
) -> tuple[int, dict[str, list[str]]]:
    dna_len = len(list(leaf_characters.values())[0])
    trees = [None] * dna_len
    full_metric = 0
    for idx in range(dna_len):
        metric, comp_tree = __small_parsimony__(
            deepcopy(tree), leaf_characters, root, idx
        )
        trees[idx] = deepcopy(comp_tree)
        full_metric += metric

    corrected_tree_names: dict[str, str] = defaultdict(str)
    for t in trees:
        for key, values in t.items():
            corrected_tree_names[key] += values
    if additional_root:
        corrected_tree_names[root] = root

    true_tree = directed2undirected(tree, corrected_tree_names)
    return full_metric, true_tree


def undirected2directed(tree: dict[str, list[str]], root: str):
    """
    Convert undirected binary tree, into directed descenent from root binary tree
    Input data:
    -----------
        tree - Undirected binary tree. Modified inplace into directed binary tree.
               Root is the first node and direction is descendant from rot to leafes.
        root - root of tree
    Output data:
    ------------
        None. Tree is modified inplace
    """
    node = root
    if is_leaf(tree, node):
        # directed_tree[node] = []
        return
    daughter, son = tree[node]
    # directed_tree[node] = [daughter, son]
    tree[son].remove(node)
    tree[daughter].remove(node)
    undirected2directed(tree, son)
    undirected2directed(tree, daughter)
    return


def directed2undirected(
    tree: dict[str, list[str]], correct_names: dict[str, str] = {}
) -> dict[str, list[str]]:
    """
    Convert directed binary tree, into undirected descenent from root binary tree
    Input data:
    -----------
        tree - Undirected binary tree.
    Output data:
    ------------
        true_tree - undirected binary tree based on tree
    """
    true_tree: dict[str, list[str]] = defaultdict(list)
    for node, adjacent_nodes in tree.items():
        for adj_node in adjacent_nodes:
            if not correct_names:
                true_tree[adj_node].append(node)
                true_tree[node].append(adj_node)
            else:
                true_tree[correct_names[adj_node]].append(correct_names[node])
                true_tree[correct_names[node]].append(correct_names[adj_node])
    return true_tree


def insert_root(tree: dict[str, list[str]], root_name: str = "root") -> str:
    try:
        # node, adj_nodes = list(tree.items())[0]
        node, adj_nodes = next(iter(tree.items()))
        adj_node = adj_nodes[0]
        tree[node].remove(adj_node)
        tree[adj_node].remove(node)
        tree[root_name] = [node, adj_node]
        tree[node].append(root_name)
        tree[adj_node].append(root_name)
        return root_name
    except Exception:
        return "Error"


def delete_root(tree: dict[str, list[str]], root_name: str = "root"):
    """
    Delete root in the undirected binary tree
    Input data:
    -----------
        tree - undrected binary tree
        root_name - key of root node in the dictionary
    Output date:
    ------------
        None. tree is modified inplace
    """
    try:
        son, daughter = tree[root_name]
        tree[son].remove(root_name)
        tree[daughter].remove(root_name)
        tree[son].append(daughter)
        tree[daughter].append(son)
        del tree[root_name]
    except KeyError:
        return "Root does not exists in the tree"
