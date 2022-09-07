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
    # tree: dict[int, dict[int, int]] = {i: {} for i in range(max_nodes)}
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

        ret = {node_i: {node_j: distance_matrix[0, 1]}, node_j: {node_i: distance_matrix[1, 0]}}
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

    #tree[max_node] = {node_i: limb_length_i, node_j: limb_length_j}
    tree[max_node][node_i] = limb_length_i
    tree[max_node][node_j] = limb_length_j
    #tree[node_i][max_node] = limb_length_i
    #tree[node_j][max_node] = limb_length_j
    tree[node_i] = {max_node: limb_length_i}
    tree[node_j] = {max_node: limb_length_j}
    return tree, max_node
