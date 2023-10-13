from typing import TextIO

import numpy as np
import numpy.typing as npt
from numba import njit


def hierarchical_clustering(
    distane_matrix: npt.NDArray[np.floating],
    f: TextIO,
) -> dict[int, dict[int, int]]:
    """
    Apply hierarchical clustering in order to get a tree, where each subcluster
    is at the same path on the tree
    Input data:
    -----------
      distane_matrix - distance matrix
    Output data:
    ------------
      Graph in the form of adjacency list
    """

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
        f.write(" ".join(map(str, map(lambda x: x + 1, new_cluster))) + "\n")
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
