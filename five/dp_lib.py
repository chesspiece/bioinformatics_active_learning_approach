import numpy as np
import numpy.typing as npt
from numba import njit

from consts.consts import BLOSUM_MATRIX, PAM_MATRIX


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


#@njit()
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


def lcs_backtrack_3d(str_1, str_2, str_3, match_reward=1, mismatch_penalty=0, indel_penalty=0):
    s = np.zeros((len(str_1) + 1, len(str_2) + 1, len(str_3) + 1), dtype=np.int64)
    backtrack = np.zeros((len(str_1) + 1, len(str_2) + 1, len(str_3) + 1), dtype=np.int64)
    for i in range(1, len(str_1) + 1):
        s[i, 0, 0] = s[i - 1, 0, 0] - indel_penalty
        backtrack[i, 0, 0] = 2
    for j in range(1, len(str_2) + 1):
        s[0, j, 0] = s[0, j - 1, 0] - indel_penalty
        backtrack[0, j, 0] = 1
    for k in range(1, len(str_3) + 1):
        s[0, 0, k] = s[0, 0, k - 1] - indel_penalty
        backtrack[0, 0, k] = 1
    for i in range(1, len(str_1) + 1):
        for j in range(1, len(str_2) + 1):
            for k in range(1, len(str_3) + 1):
                #match1 = match_reward if str_1[i - 1] == str_2[j - 1] else -mismatch_penalty
                #match2 = match_reward if str_1[i - 1] == str_3[j - 1] else -mismatch_penalty
                #match3 = match_reward if str_2[i - 1] == str_3[j - 1] else -mismatch_penalty

                match4 = match_reward if str_1[i - 1] == str_2[j - 1] and \
                         str_2[j - 1] == str_3[k - 1] else -mismatch_penalty

                s[i, j] = max(
                    s[i - 1, j, k] - indel_penalty,
                    s[i, j - 1, k] - indel_penalty,
                    s[i, j, k - 1] - indel_penalty,
                    s[i - 1, j - 1, k] - indel_penalty,
                    s[i - 1, j, k - 1] - indel_penalty,
                    s[i, j - 1, k - 1] - indel_penalty,
                    s[i - 1, j - 1, k - 1] + match4
                )
                """
                if s[i, j] == s[i - 1, j] - indel_penalty:
                    backtrack[i, j] = 2
                elif s[i, j] == s[i, j - 1] - indel_penalty:
                    backtrack[i, j] = 1
                elif match == match_reward:
                    backtrack[i, j] = 3
                else:
                    backtrack[i, j] = 4
                """
    return s[len(str_1), len(str_2), len(str_3)] #backtrack, s[len(str_1), len(str_2)]



@njit()
def middle_node(str_1, str_2, match_reward=1, mismatch_penalty=0, indel_penalty=0):
    middle = len(str_2) // 2
    reversed_str_1 = str_1[::-1]

    reversed_str_2_middle = str_2[middle:][::-1]
    str_2_middle = str_2[:middle]

    def comp_i_col(str_1, str_2, match_reward, mismatch_penalty, indel_penalty):
        score = np.zeros((len(str_1) + 1, 2), dtype=np.int64)
        for i in range(1, len(str_1) + 1):
            score[i, 0] = score[i - 1, 0] - indel_penalty
        #score[0, 1] = score[0, 0] - indel_penalty

        for j in range(1, len(str_2) + 1):
            score[0, j % 2] = score[0, (j + 1) % 2] - indel_penalty
            for i in range(1, len(str_1) + 1):
                match = (
                    match_reward if str_1[i - 1] == str_2[j - 1] else -mismatch_penalty
                )
                score[i, j % 2] = max(
                    score[i - 1, j % 2] - indel_penalty,
                    score[i, (j - 1) % 2] - indel_penalty,
                    score[i - 1, (j - 1) % 2] + match,
                )
        if len(str_2) % 2 == 0:
            score = score[:, ::-1]

        return score

    so_m = comp_i_col(
        str_1, str_2_middle, match_reward, mismatch_penalty, indel_penalty
    )
    m_si = np.flip(
        comp_i_col(
            reversed_str_1,
            reversed_str_2_middle,
            match_reward,
            mismatch_penalty,
            indel_penalty,
        )
    )
    edge1 = (np.argmax(so_m[:, -1] + m_si[:, 0]), middle)
    f_node_i, f_node_j = edge1

    #print((str_1))
    #print((str_2))
    #print(len(str_1))
    #print(len(str_2))
    #print(f_node_i)
    #print(m_si)
    #print(f_node_i)
    #print(m_si)
    sz1, _ = m_si.shape
    #print(sz1)
    #print(f_node_i)
    #if f_node_i == sz1 - 1:
    ##    s_node_i = f_node_i
    #    s_node_j = f_node_j + 1
    #    edge_dir = 1
    #    return (int(f_node_i), f_node_j, int(s_node_i), s_node_j), edge_dir
    if f_node_i < sz1 - 1:
        if m_si[f_node_i, 0] + indel_penalty == m_si[f_node_i, 1]:
            s_node_i = f_node_i
            s_node_j = f_node_j + 1
            edge_dir = 1
        elif m_si[f_node_i, 0] + indel_penalty == m_si[f_node_i + 1, 0]:
            s_node_i = f_node_i + 1
            s_node_j = f_node_j
            edge_dir = 2
        else:
            s_node_i = f_node_i + 1
            s_node_j = f_node_j + 1
            edge_dir = 3
    else:
            s_node_i = f_node_i
            s_node_j = f_node_j + 1
            edge_dir = 1
    return (int(f_node_i), f_node_j, int(s_node_i), s_node_j), edge_dir


def linear_space_aligment(
    str_1: str,
    str_2: str,
    top: int,
    bottom: int,
    left: int,
    right: int,
    path: list[int],
    math_reward: int = 1,
    mismatch_penalty: int = 0,
    indel_penalty: int = 0,
):
    if left == right:
        new_path = [2] * (bottom - top)
        path += new_path
        return
    if top == bottom:
        new_path = [1] * (right - left)
        path += new_path
        return
    middle = int(np.floor((right + left) / 2))
    print("Yay1")
    middle_edge_coord, edge_direction = middle_node(
        str_1[top:bottom], str_2[left:right], math_reward, mismatch_penalty, indel_penalty
    )
    mid_node, _, _, _ = middle_edge_coord
    mid_node += top#left
    print("Yay2")
    linear_space_aligment(
        str_1,
        str_2,
        top,
        mid_node,
        left,
        middle,
        path,
        math_reward,
        mismatch_penalty,
        indel_penalty,
    )
    path.append(edge_direction)
    if edge_direction == 6:
        path += [1]*(right - middle)
        return
    if edge_direction == 1 or edge_direction == 3:
        middle += 1
    if edge_direction == 2 or edge_direction == 3:
        mid_node += 1
    linear_space_aligment(
        str_1,
        str_2,
        mid_node,
        bottom,
        middle,
        right,
        path,
        math_reward,
        mismatch_penalty,
        indel_penalty,
    )
    return


def find_path_back(dna_1: str, dna_2: str, path: list[int], match_reward=1, mismatch_penalty=0, indel_penalty=0):
    i = len(dna_1) - 1
    j = len(dna_2) - 1
    aligment_1: str = ""
    aligment_2: str = ""
    score = 0
    for edge_dir in reversed(path):
        if edge_dir == 3:
            score += match_reward if dna_1[i] == dna_2[j] else -mismatch_penalty
            aligment_1 += dna_1[i] 
            aligment_2 += dna_2[j]
            i -= 1
            j -= 1
        if edge_dir == 2:
            score -= indel_penalty
            aligment_1 += dna_1[i] 
            aligment_2 += "-"#dna_2[j]
            i -= 1
            #j -= 1
        if edge_dir == 1:
            score -= indel_penalty
            aligment_1 += "-"#dna_1[i] 
            aligment_2 += dna_2[j]
            #i -= 1
            j -= 1
    return score, aligment_1[::-1], aligment_2[::-1]




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
            if s[i, j] > max_score:
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
        s[i, 0] = 0  # s[i - 1, 0] - indel_penalty
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
            if s[i, j] > max_score:
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


def lcs_backtrack_overlap(
    str_1, str_2, match_reward=1, mismatch_penalty=1, indel_penalty=1
):
    s = np.zeros((len(str_1) + 1, len(str_2) + 1), dtype=np.int64)
    backtrack = np.zeros((len(str_1) + 1, len(str_2) + 1), dtype=np.int64)
    max_score = -np.inf
    max_i, max_j = 0, 0
    for i in range(1, len(str_1) + 1):
        s[i, 0] = 0  # s[i - 1, 0] - indel_penalty
        backtrack[i, 0] = 0
    for i in range(1, len(str_2) + 1):
        s[0, i] = 0  # s[0, i - 1] - indel_penalty
        backtrack[0, i] = 0
    for i in range(1, len(str_1) + 1):
        for j in range(1, len(str_2) + 1):
            match = match_reward if str_1[i - 1] == str_2[j - 1] else -mismatch_penalty
            s[i, j] = max(
                s[i - 1, j] - indel_penalty,
                s[i, j - 1] - indel_penalty,
                s[i - 1, j - 1] + match,
            )
            if s[i, j] > max_score:
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
    # backtrack, score = lcs_backtrack(
    #   str1, str2, match_reward, mismatch_penalty, indel_penalty
    # )
    backtrack, score, max_i, max_j = lcs_backtrack_overlap(
        str1, str2, match_reward, mismatch_penalty, indel_penalty
    )
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
        # __output_lcs1(backtrack, str1, len(str1), len(str2)),
        # __output_lcs2(backtrack, str2, len(str1), len(str2)),
        # __output_lcs0(backtrack, str2, len(str1), len(str2)),
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
