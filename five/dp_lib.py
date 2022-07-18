import numpy as np
import numpy.typing as npt
from numba import njit


def dp_change(money: int, coins: npt.NDArray[np.int64]) -> int:
    min_num_coins = np.zeros(money + 1, dtype=np.int64)
    min_num_coins[0] = 0
    for money_sum in range(1, min_num_coins.size):
        min_num_coins[money_sum] = np.iinfo(np.int64).max
        for coin in coins:
            if money_sum >= coin:
                candidate = min_num_coins[money_sum - coin] + 1
                if candidate < min_num_coins[money_sum]:
                    min_num_coins[money_sum] = candidate
    return min_num_coins[money]


@njit()
def manhattan_tourist(n: int, m: int, down, right) -> int:
    s = np.zeros((n, m), dtype=np.int64)
    for i in range(1, n):
        s[i, 0] = s[i - 1, 0] + down[i - 1, 0]
        s[0, i] = s[0, i - 1] + right[0, i - 1]
    for i in range(1, n):
        for j in range(1, m):
            s[i, j] = max(s[i - 1, j] + down[i - 1, j], s[i, j - 1] + right[i, j - 1])
    return s[n - 1, m - 1]


def lcs_backtrack(str_1, str_2, match_reward=1, mismatch_penalty=0, indel_penalty=0):
    s = np.zeros((len(str_1), len(str_2)), dtype=np.int64)
    backtrack = np.zeros((len(str_1), len(str_2)), dtype=np.int64)
    for i in range(len(str_1)):
        s[i, 0] = 0
    for i in range(len(str_2)):
        s[0, i] = 0
    for i in range(1, len(str_1)):
        for j in range(1, len(str_2)):
            match = 1 if str_1[i] == str_2[j] else 0
            s[i, j] = max(s[i - 1, j], s[i, j - 1], s[i - 1, j - 1] + match)
            if s[i, j] == s[i - 1, j]:
                backtrack[i, j] = 1
            elif s[i, j] == s[i, j - 1]:
                backtrack[i, j] = 2
            else:
                backtrack[i, j] = 3
    return backtrack


def output_lcs(str1, str2):
    backtrack = lcs_backtrack(str1, str2)
    print(backtrack)

    def __output_lcs(backtrack, str_1, i, j):
        if i == -1 or j == -1:
            return ""
        if backtrack[i, j] == 1:
            return __output_lcs(backtrack, str_1, i - 1, j)
        elif backtrack[i, j] == 2:
            return __output_lcs(backtrack, str_1, i, j - 1)
        else:
            return __output_lcs(backtrack, str_1, i - 1, j - 1) + str_1[i]

    return __output_lcs(backtrack, str1, len(str1) - 1, len(str2) - 1)


class Grph_path:
    def __init__(self, V_quant, E_quant):
        self.V = V_quant
        self.E = E_quant
        self.edges_list = [[]] * (self.V)
        self.reverse_edges_list = [[]] * (self.V)

    def build_graph(self, edges_list):
        for i in range(0, self.E):
            vertex1, vertex2, weight = edges_list[i]
            self.edges_list[vertex1] = self.edges_list[vertex1] + [(vertex2, weight)]
            self.reverse_edges_list[vertex2] = self.reverse_edges_list[vertex2] + [
                (vertex1, weight)
            ]

    def topologically_sort(self):
        """
        Not actually topologicall sort.
        Graph has cycles in this task
        """
        sorted_grph = []
        visited_list = [False] * self.V
        for v_comp in range(0, self.V):
            if visited_list[v_comp] is True:
                continue
            vertex = v_comp
            stack_l = []
            stack_l.append(vertex)
            while len(stack_l) > 0:
                vrtx = stack_l[len(stack_l) - 1]
                if visited_list[vrtx] is False:
                    visited_list[vrtx] = True
                else:
                    stack_l.pop()
                    if not vrtx in sorted_grph:
                        sorted_grph.append(vrtx)
                    continue
                for i, _ in self.edges_list[vrtx]:
                    if visited_list[i] is False:
                        stack_l.append(i)
        return sorted_grph[::-1]

    def longest_path(self, source, sink):
        """
        TODO: source not can only have indegree=0 now!!
        """
        path_length = [-np.inf] * self.V
        path_pred = [-1] * self.V
        path_length[source] = 0
        top_grph = self.topologically_sort()
        for vertex in top_grph:
            for vertex_source, weight in self.reverse_edges_list[vertex]:
                if path_length[vertex] < path_length[vertex_source] + weight:
                    path_length[vertex] = path_length[vertex_source] + weight
                    path_pred[vertex] = vertex_source
        path = []
        v = sink
        while v != source:
            path.append(v)
            v = path_pred[v]
        path.append(v)
        return path_length[sink], reversed(path)
