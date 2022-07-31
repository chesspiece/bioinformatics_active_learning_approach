from collections import defaultdict, deque

import numpy as np
import numpy.typing as npt
from numba import njit, prange
from numba.typed import List


@njit()
def greedy_sorting(perm: npt.NDArray[np.int64]) -> tuple[int, list[list[int]]]:
    approx_dist: int = 0
    reversal_sequence: list[list[int]] = []
    for idx, perm_elem in enumerate(perm):
        if np.abs(perm_elem) != idx + 1:
            elem_pos = np.where(np.abs(perm) == (idx + 1))[0][0]
            perm[idx : elem_pos + 1] = -np.flip(perm[idx : elem_pos + 1])
            reversal_sequence.append(list(perm))
            approx_dist += 1
        if perm[idx] == -(idx + 1):
            perm[idx] = -perm[idx]
            reversal_sequence.append(list(perm))
            approx_dist += 1
    return approx_dist, reversal_sequence


def breakpoints_number(perm: npt.NDArray[np.int64]) -> int:
    breakpoints_quant: int = 0 if perm[0] == 1 else 1
    for idx in range(1, len(perm)):
        if perm[idx] - 1 != perm[idx - 1]:
            breakpoints_quant += 1
    if perm[-1] != len(perm):
        breakpoints_quant += 1
    return breakpoints_quant


#@njit(parallel=True)
def chromosome_to_cycle(chromosome: list[int]) -> list[int]:
    nodes = [0] * (2 * len(chromosome))
    for idx in prange(len(chromosome)):
        synteny = chromosome[idx]
        if synteny > 0:
            nodes[2 * (idx + 1) - 1 - 1] = 2 * synteny - 1
            nodes[2 * (idx + 1) - 1] = 2 * synteny
        else:
            nodes[2 * (idx + 1) - 1 - 1] = -2 * synteny
            nodes[2 * (idx + 1) - 1] = -2 * synteny - 1
    return nodes


#@njit(parallel=True)
def cycle_to_chromosome(nodes: list[int]) -> list[int]:
    chromosome = [0] * (len(nodes) // 2)
    for idx in prange(len(nodes) // 2):
        if nodes[2 * (idx + 1) - 1 - 1] < nodes[2 * (idx + 1) - 1]:
            chromosome[idx] = nodes[2 * (idx + 1) - 1] // 2
        else:
            chromosome[idx] = -nodes[2 * (idx + 1) - 1] // 2
    return chromosome


def colored_edges(chromosome_list: list[list[int]]) -> list[tuple[int, int]]:
    edges: list[tuple[int, int]] = []
    for chromosome in chromosome_list:
        nodes = chromosome_to_cycle(List(chromosome))
        for idx in range(len(chromosome) - 1):
            edges.append((nodes[2 * idx + 1], nodes[2 * idx + 2]))
        edges.append((nodes[-1], nodes[0]))
    return edges


def different_colored_edges(
    chromosome_list: list[list[int]], color: int = 0
) -> list[tuple[int, int, int]]:
    edges: list[tuple[int, int, int]] = []
    for chromosome in chromosome_list:
        nodes = chromosome_to_cycle(List(chromosome))
        for idx in range(len(chromosome) - 1):
            edges.append((nodes[2 * idx + 1], nodes[2 * idx + 2], color))
        edges.append((nodes[-1], nodes[0], color))
    return edges


def ordered_edges_to_cycle(edges: list[tuple[int, int]]) -> list[list[int]]:
    cycles = []
    curr_cycle = []
    for i, j in edges:
        if j < i:
            cycles.append([j] + curr_cycle + [i])
            curr_cycle = []
        else:
            curr_cycle += [i, j]
    return cycles


def edges_to_genome(edges: list[tuple[int, int]]) -> list[list[int]]:
    chromosomes: list[list[int]] = []
    black_edges = uncolored_2_colored(generate_black_edges(len(edges)), color=0)
    colored_edges = uncolored_2_colored(edges, color=1)
    adj_list = create_adjacency_list(black_edges + colored_edges)
    cycles = dfs_find_cicle_black(adj_list)
    for cycle in cycles:
        chromosomes.append(cycle_to_chromosome(List(cycle)))
    return chromosomes


def two_breaks_genome_graph(
    genome_edges: list[tuple[int, int]], i: int, j: int, k: int, l: int
) -> list[tuple[int, int]]:
    if (i, j) in genome_edges:
        genome_edges.remove((i, j))
    else:
        genome_edges.remove((j, i))
    if (k, l) in genome_edges:
        genome_edges.remove((k, l))
    else:
        genome_edges.remove((l, k))
    genome_edges += [(i, k), (j, l)]
    genome_edges = sorted(genome_edges, key=lambda x: x[0])
    return genome_edges


def two_breaks_genome(
    genome: list[list[int]], i: int, j: int, k: int, l: int
) -> list[list[int]]:
    genome_graph = colored_edges(genome)
    genome_graph = two_breaks_genome_graph(genome_graph, i, j, k, l)
    return edges_to_genome(genome_graph)


def generate_black_edges(num_syn: int) -> list[tuple[int, int]]:
    return [(i, i + 1) for i in range(1, 2 * num_syn, 2)]


def uncolored_2_colored(
    edges: list[tuple[int, int]], color: int = 0
) -> list[tuple[int, int, int]]:
    return [(x, y, color) for x, y in edges]


def colored_2_uncolored(edges: list[tuple[int, int, int]]) -> list[tuple[int, int]]:
    return [(x, y) for x, y, _ in edges]


def create_adjacency_list(
    chromosome_list: list[tuple[int, int, int]]
) -> dict[int, list[tuple[int, int]]]:
    adj_list: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for i, j, color in chromosome_list:
        adj_list[i] += [(j, color)]
        adj_list[j] += [(i, color)]
    return adj_list


def dfs_find_cicle_black(adjacency_list: dict[int, list[tuple[int, int]]]):
    visited_node = defaultdict(bool)
    circles_list = []
    for curr_node in sorted(adjacency_list.keys()):
        if visited_node[curr_node] or (curr_node % 2 != 1):
            continue
        v1, v2 = adjacency_list[curr_node]
        if v1[0] == curr_node + 2 or v1[0] == curr_node + 3:
            curr_node = v2[0]
        elif v2[0] == curr_node + 2 or v2[0] == curr_node + 3:
            curr_node = v1[0]
        vertexes = [curr_node]
        circle = []
        curr_color = 1
        while vertexes:
            curr_node = vertexes.pop()
            circle.append(curr_node)
            visited_node[curr_node] = True
            for node, color in adjacency_list[curr_node]:
                if not visited_node[node] and color != curr_color:
                    vertexes.append(node)
                elif color != curr_color:
                    # circle.append(node)
                    circles_list.append(circle)
                    break
            curr_color = 0 if curr_color == 1 else 1
    return circles_list


def dfs_find_cicle_colored(adjacency_list: dict[int, list[tuple[int, int]]]):
    visited_node = defaultdict(bool)
    circles_list = []
    for curr_node in sorted(adjacency_list.keys()):
        if visited_node[curr_node]:
            continue
        vertexes = [curr_node]
        circle = []
        curr_color = 1
        while vertexes:
            curr_node = vertexes.pop()
            circle.append(curr_node)
            visited_node[curr_node] = True
            for node, color in adjacency_list[curr_node]:
                if not visited_node[node] and color != curr_color:
                    vertexes.append(node)
                elif color != curr_color:
                    # circle.append(node)
                    circles_list.append(circle)
                    break
            curr_color = 0 if curr_color == 1 else 1
    return circles_list


def shortest_rearrangmen(
    genome_1: list[list[int]], genome_2: list[list[int]]
) -> list[str]:
    tmp_1 = [" ".join([f"{gen:+d}" for gen in chromosome]) for chromosome in genome_1]
    tmp = "(" + ")(".join(tmp_1) + ")"
    output_rear = [tmp]
    red_edges = colored_edges(genome_1)
    blue_edges = colored_edges(genome_2)
    
    red_edges = different_colored_edges(genome_1, color=1)
    blue_edges = different_colored_edges(genome_2, color=0)
 
    #cycles = dfs_find_cicle_black(create_adjacency_list(red_edges + blue_edges))
    cycles = dfs_find_cicle_colored(create_adjacency_list(red_edges + blue_edges))

    for cycle in cycles:
        if len(cycles) == 2:
            continue
        for i in range(2, len(cycle), 2):
            break_points = [cycle[i - 2], cycle[-1], cycle[i - 1], cycle[i]]
            genome_1 = two_breaks_genome(genome_1, *break_points)
            tmp_1 = [" ".join([f"{gen:+d}" for gen in chromosome]) for chromosome in genome_1]
            tmp = "(" + ")(".join(tmp_1) + ")"
            output_rear.append(tmp)
    return output_rear
