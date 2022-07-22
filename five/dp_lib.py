import numpy as np
import numpy.typing as npt
from numba import njit

from consts.consts import PAM_MATRIX, BLOSUM_MATRIX


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
    s = np.zeros((len(str_1) + 1, len(str_2) + 1), dtype=np.int64)
    backtrack = np.zeros((len(str_1) + 1, len(str_2) + 1), dtype=np.int64)
    for i in range(1, len(str_1) + 1):
        s[i, 0] = s[i - 1, 0] - indel_penalty
        backtrack[i, 0] = 2
    for i in range(1, len(str_2) + 1):
        s[0, i] = s[0, i - 1] - indel_penalty
        backtrack[0, i] = 1
    for i in range(1, len(str_1) + 1):
        for j in range(1, len(str_2) + 1):
            match = match_reward if str_1[i - 1] == str_2[j - 1] else -mismatch_penalty
            s[i, j] = max(
                s[i - 1, j] - indel_penalty,
                s[i, j - 1] - indel_penalty,
                s[i - 1, j - 1] + match,
            )
            if s[i, j] == s[i - 1, j] - indel_penalty:
                backtrack[i, j] = 2
            elif s[i, j] == s[i, j - 1] - indel_penalty:
                backtrack[i, j] = 1
            elif match == match_reward:
                backtrack[i, j] = 3
            else:
                backtrack[i, j] = 4
    return backtrack, s[len(str_1), len(str_2)]


def lcs_backtrack_local(str_1, str_2, score_matrix=PAM_MATRIX, indel_penalty=5):
    s = np.zeros((len(str_1) + 1, len(str_2) + 1), dtype=np.int64)
    backtrack = np.zeros((len(str_1) + 1, len(str_2) + 1), dtype=np.int64)
    max_score = -np.inf
    max_i, max_j = 0, 0
    for i in range(1, len(str_1) + 1):
        s[i, 0] = s[i - 1, 0] - indel_penalty
        backtrack[i, 0] = 2
    for i in range(1, len(str_2) + 1):
        s[0, i] = s[0, i - 1] - indel_penalty
        backtrack[0, i] = 1
    for i in range(1, len(str_1) + 1):
        for j in range(1, len(str_2) + 1):
            # match = match_reward if str_1[i - 1] == str_2[j - 1] else -mismatch_penalty
            match = score_matrix[str_1[i - 1]][str_2[j - 1]]
            s[i, j] = max(
                0,
                s[i - 1, j] - indel_penalty,
                s[i, j - 1] - indel_penalty,
                s[i - 1, j - 1] + match,
            )
            if s[i,  j] > max_score:
                max_score = s[i, j]
                max_i, max_j = i, j
            if s[i, j] == s[i - 1, j] - indel_penalty:
                backtrack[i, j] = 2
            elif s[i, j] == s[i, j - 1] - indel_penalty:
                backtrack[i, j] = 1
            elif (s[i, j] == s[i - 1, j - 1] + match) and str_1[i - 1] == str_2[j - 1]:
                backtrack[i, j] = 3
            elif (s[i, j] == s[i - 1, j - 1] + match) and str_1[i - 1] != str_2[j - 1]:
                backtrack[i, j] = 4
            else:
                backtrack[i, j] = 0
    return backtrack, s[max_i, max_j], max_i, max_j


def lcs_backtrack_fitting(str_1, str_2, score_matrix=BLOSUM_MATRIX, indel_penalty=1):
    s = np.zeros((len(str_1) + 1, len(str_2) + 1), dtype=np.int64)
    backtrack = np.zeros((len(str_1) + 1, len(str_2) + 1), dtype=np.int64)
    max_score = -np.inf
    max_i, max_j = 0, 0
    for i in range(1, len(str_1) + 1):
        s[i, 0] = 0# s[i - 1, 0] - indel_penalty
        backtrack[i, 0] = 0
    for i in range(1, len(str_2) + 1):
        s[0, i] = s[0, i - 1] - indel_penalty
        backtrack[0, i] = 1
    for i in range(1, len(str_1) + 1):
        for j in range(1, len(str_2) + 1):
            # match = match_reward if str_1[i - 1] == str_2[j - 1] else -mismatch_penalty
            match = score_matrix[str_1[i - 1]][str_2[j - 1]]
            s[i, j] = max(
                s[i - 1, j] - indel_penalty,
                s[i, j - 1] - indel_penalty,
                s[i - 1, j - 1] + match,
            )
            if s[i,  j] > max_score:
                max_score = s[i, j]
                max_i, max_j = i, j
            if s[i, j] == s[i - 1, j] - indel_penalty:
                backtrack[i, j] = 2
            elif s[i, j] == s[i, j - 1] - indel_penalty:
                backtrack[i, j] = 1
            elif (s[i, j] == s[i - 1, j - 1] + match) and str_1[i - 1] == str_2[j - 1]:
                backtrack[i, j] = 3
            elif (s[i, j] == s[i - 1, j - 1] + match) and str_1[i - 1] != str_2[j - 1]:
                backtrack[i, j] = 4
            else:
                backtrack[i, j] = 0
    print(s)
    return backtrack, s[max_i, max_j], max_i, max_j


def lcs_backtrack_overlap(str_1, str_2, match_reward=1,  mismatch_penalty=1, indel_penalty=1):
    s = np.zeros((len(str_1) + 1, len(str_2) + 1), dtype=np.int64)
    backtrack = np.zeros((len(str_1) + 1, len(str_2) + 1), dtype=np.int64)
    max_score = -np.inf
    max_i, max_j = 0, 0
    for i in range(1, len(str_1) + 1):
        s[i, 0] = 0# s[i - 1, 0] - indel_penalty
        backtrack[i, 0] = 0
    for i in range(1, len(str_2) + 1):
        s[0, i] = 0 #s[0, i - 1] - indel_penalty
        backtrack[0, i] = 0
    for i in range(1, len(str_1) + 1):
        for j in range(1, len(str_2) + 1):
            match = match_reward if str_1[i - 1] == str_2[j - 1] else -mismatch_penalty
            s[i, j] = max(
                s[i - 1, j] - indel_penalty,
                s[i, j - 1] - indel_penalty,
                s[i - 1, j - 1] + match,
            )
            if s[i,  j] > max_score:
                max_score = s[i, j]
                max_i, max_j = i, j
            if s[i, j] == s[i - 1, j] - indel_penalty:
                backtrack[i, j] = 2
            elif s[i, j] == s[i, j - 1] - indel_penalty:
                backtrack[i, j] = 1
            elif (s[i, j] == s[i - 1, j - 1] + match) and str_1[i - 1] == str_2[j - 1]:
                backtrack[i, j] = 3
            elif (s[i, j] == s[i - 1, j - 1] + match) and str_1[i - 1] != str_2[j - 1]:
                backtrack[i, j] = 4
            else:
                backtrack[i, j] = 0
    print(s)
    return backtrack, s[max_i, max_j], max_i, max_j



def output_lcs(str1, str2, match_reward=1, mismatch_penalty=0, indel_penalty=0):
    #backtrack, score = lcs_backtrack(
    #   str1, str2, match_reward, mismatch_penalty, indel_penalty
    #)
    backtrack, score, max_i, max_j = lcs_backtrack_overlap(str1, str2, match_reward, mismatch_penalty, indel_penalty)
    print(backtrack)

    def __output_lcs0(backtrack, str_1, i, j):
        if i == 0 and j == 0:
            return ""
        if backtrack[i, j] == 0:
            return ""
        elif backtrack[i, j] == 1:
            return __output_lcs0(backtrack, str_1, i, j - 1)
        elif backtrack[i, j] == 2:
            return __output_lcs0(backtrack, str_1, i - 1, j)
        elif backtrack[i, j] == 3:
            return __output_lcs0(backtrack, str_1, i - 1, j - 1) + str_1[i - 1]
        else:
            return __output_lcs0(backtrack, str_1, i - 1, j - 1)

    def __output_lcs1(backtrack, str_1, i, j):
        if i == 0 and j == 0:
            return ""
        if backtrack[i, j] == 0:
            return ""
        elif backtrack[i, j] == 1:
            return __output_lcs1(backtrack, str_1, i, j - 1) + "-"
        elif backtrack[i, j] == 2:
            return __output_lcs1(backtrack, str_1, i - 1, j) + str_1[i - 1]
        else:
            return __output_lcs1(backtrack, str_1, i - 1, j - 1) + str_1[i - 1]

    def __output_lcs2(backtrack, str_1, i, j):
        if j == 0 and i == 0:
            return ""
        if backtrack[i, j] == 0:
            return ""
        elif backtrack[i, j] == 1:
            return __output_lcs2(backtrack, str_1, i, j - 1) + str_1[j - 1]
        elif backtrack[i, j] == 2:
            return __output_lcs2(backtrack, str_1, i - 1, j) + "-"
        else:
            return __output_lcs2(backtrack, str_1, i - 1, j - 1) + str_1[j - 1]

    return (
        score,
        __output_lcs1(backtrack, str1, max_i, max_j),
        __output_lcs2(backtrack, str2, max_i, max_j),
        #__output_lcs1(backtrack, str1, len(str1), len(str2)),
        #__output_lcs2(backtrack, str2, len(str1), len(str2)),
        #__output_lcs0(backtrack, str2, len(str1), len(str2)),
    )


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
