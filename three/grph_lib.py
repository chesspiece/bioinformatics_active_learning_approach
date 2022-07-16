from collections import defaultdict
from copy import deepcopy
from typing import Tuple

from numba import njit


class myGrph:
    def __init__(self, vertex_list: list[str]):
        self.vertex_list = vertex_list
        self.V_quant = len(vertex_list)
        self.edges_list: dict[str, list[str]] = defaultdict(list)
        self.__fill_connections__()

    def __fill_connections__(self):
        for i in range(self.V_quant):
            for j in range(self.V_quant):
                if i == j:
                    continue
                if self.vertex_list[i][1::] == self.vertex_list[j][:-1]:
                    self.edges_list[self.vertex_list[i]] += [self.vertex_list[j]]

    def write_edges(self, filename: str):
        with open(filename, "w") as f:
            for key, value in self.edges_list.items():
                f.write(key + ":" + " " + " ".join(value) + "\n")


# @njit()
def str_comp(dna: str, k: int) -> list[str]:
    """
    decompose string into k-mers
    """
    res = []
    for i in range(0, len(dna) - k + 1):
        res.append(dna[i : i + k])
    return res


def write_edges(edges_list: dict[str, list[str]], filename: str) -> None:
    with open(filename, "w") as f:
        for key, value in edges_list.items():
            f.write(key + ":" + " " + " ".join(value) + "\n")


def find_euler_cycle(edges_list, v):
    path = [v]
    while edges_list[v]:
        path = path[:1] + find_euler_cycle(edges_list, edges_list[v].pop()) + path[1:]
    return path


def find_euler_path(edges_list, i_deg, o_deg, v_list):
    start_node = find_initial_node(v_list, i_deg, o_deg)
    return find_euler_cycle(edges_list, start_node)


def find_initial_node(
    v_list: list[str],
    i_deg: dict[str, int],
    o_deg: dict[str, int],
):
    # initial_node = None
    initial_node = v_list[0]
    for node in v_list:
        if i_deg[node] < o_deg[node]:
            initial_node = node
            break
    return initial_node


def max_nonbranch_paths(e_list, v_list, i_deg, o_deg):
    paths = []
    checked = []
    for v in v_list:
        if not (i_deg[v] == 1 and o_deg[v] == 1):
            if o_deg[v] > 0:
                checked.append(v)
                for u in e_list[v]:
                    path = [v, u]
                    while i_deg[u] == 1 and o_deg[u] == 1:
                        checked.append(v)
                        path.append(e_list[u][0])
                        u = e_list[u][0]
                    paths.append(path)
    for v in v_list:
        if (i_deg[v] == 1 and o_deg[v] == 1 and (not v in checked)):
            for u in e_list[v]:
                path = [v, u]
                checked.append(v)
                while i_deg[u] == 1 and o_deg[u] == 1 and u != v:
                    checked.append(u)
                    path.append(e_list[u][0])
                    u = e_list[u][0]
                if u == v:
                    paths.append(path)
    return paths


def helper_dna_path(pth):
    dna = pth[0]
    dna_append = []
    for k_mer in pth[1::]:
        dna_append.append(k_mer[-1])
    return dna + "".join(dna_append)


class deBrujin_grph:
    def __init__(self, dna_list: list[str]):
        self.vertex_list = []
        self.V_quant = 0
        self.edges_list: dict[str, list[str]] = defaultdict(list)
        self.i_deg: dict[str, int] = defaultdict(int)
        self.o_deg: dict[str, int] = defaultdict(int)
        self.__fill_connections__(dna_list)

    def __fill_connections__(self, dna_list: list[str]):
        for i in range(len(dna_list)):
            self.i_deg[dna_list[i][1::]] += 1
            self.o_deg[dna_list[i][:-1]] += 1
            if not (dna_list[i][1::] in self.vertex_list):
                self.vertex_list.append(dna_list[i][1::])
            if not (dna_list[i][:-1] in self.vertex_list):
                self.vertex_list.append(dna_list[i][:-1])
            self.edges_list[dna_list[i][:-1]] += [dna_list[i][1::]]
        self.V_quant = len(self.vertex_list)

    def write_edges(self, filename: str):
        with open(filename, "w") as f:
            for key, value in self.edges_list.items():
                f.write(key + ":" + " " + " ".join(value) + "\n")

    def find_dna(self):
        answ = find_euler_path(
            deepcopy(self.edges_list), self.i_deg, self.o_deg, self.vertex_list
        )
        dna = answ[0]
        dna_append = []
        for k_mer in answ[1::]:
            dna_append.append(k_mer[-1])
        return dna + "".join(dna_append)

    def find_dna_cycle(self):
        answ = find_euler_cycle(deepcopy(self.edges_list), self.vertex_list[0])
        dna = answ[0]
        k = len(dna)
        dna_append = []
        for k_mer in answ[1::]:
            dna_append.append(k_mer[-1])
        ret = dna + "".join(dna_append)
        return ret[0 : len(ret) - k]


class deBrujin_grph2:
    def __init__(self, dna_list: list[str]):
        self.vertex_list = []
        self.V_quant = 0
        self.edges_list: dict[tuple[str, str], list[tuple[str, str]]] = defaultdict(
            list
        )
        self.i_deg: dict[tuple[str, str], int] = defaultdict(int)
        self.o_deg: dict[tuple[str, str], int] = defaultdict(int)
        self.__fill_connections__(dna_list)

    def __fill_connections__(self, dna_list: list[str]):
        for i in range(len(dna_list)):
            vrtx = dna_list[i].split("|")
            vrtx1 = vrtx[0]
            vrtx2 = vrtx[1]

            self.i_deg[(vrtx1[1::], vrtx2[1::])] += 1
            self.o_deg[(vrtx1[:-1], vrtx2[:-1])] += 1

            if not (dna_list[i][1::] in self.vertex_list):
                self.vertex_list.append((vrtx1[1::], vrtx2[1::]))

            if not (dna_list[i][:-1] in self.vertex_list):
                self.vertex_list.append((vrtx1[:-1], vrtx2[:-1]))

            self.edges_list[(vrtx1[:-1], vrtx2[:-1])] += [(vrtx1[1::], vrtx2[1::])]
        self.V_quant = len(self.vertex_list)

    def find_dna(self, d):
        answ = find_euler_path(
            deepcopy(self.edges_list), self.i_deg, self.o_deg, self.vertex_list
        )
        dna1, dna2 = answ[0]
        k = len(dna1)
        dna_append1 = []
        dna_append2 = [dna2]
        for k_mer in answ[1::]:
            t1, t2 = k_mer
            dna_append1.append(t1[-1])
            dna_append2.append(t2[-1])
        return (
            dna1
            + "".join(dna_append1)
            + "".join(dna_append2[len(dna_append2) - d - k - 1 : len(dna_append2)])
        )

    def find_dna_cycle(self):
        answ = find_euler_cycle(deepcopy(self.edges_list), self.vertex_list[0])
        dna, _ = answ[0]
        # k = len(dna)
        dna_append = []
        for k_mer in answ[1::]:
            t1, _ = k_mer
            dna_append.append(t1[-1])
        ret = dna + "".join(dna_append)
        # return ret[0 : len(ret) - k]
        return ret


class Euler_grph:
    def __init__(self, input_filename: str):
        self.vertex_list: list[str] = []
        self.V_quant = 0
        self.edges_list: dict[str, list[str]] = defaultdict(list)
        self.i_deg: dict[str, int] = defaultdict(int)
        self.o_deg: dict[str, int] = defaultdict(int)
        self.__fill_connections__(input_filename)

    def __fill_connections__(self, input_filename):
        with open(input_filename, "r") as f:
            for line in f:
                edges = line.strip().split(":")
                self.vertex_list.append(edges[0])
                self.V_quant += 1
                out_edges = edges[1].strip().split()
                self.o_deg[edges[0]] += len(out_edges)
                for edg in out_edges:
                    self.i_deg[edg] += 1
                self.edges_list[edges[0]] += out_edges
