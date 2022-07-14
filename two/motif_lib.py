from copy import deepcopy

import numpy as np
import numpy.typing as npt
from numba import njit, prange
from numba.typed import List

from one.cs_excercise import neighbours
from one.skew import hamming_str


@njit()
def motif_enum(dna: list[str], k: int, d: int) -> list[str]:
    """
    Find if motif presented in all dna strings or substrings with at most 2d mismatches
    Very slow implementation
    Input data:
    -----------
        dna - input dna string
        d - allowed number of mismathces
        k - k-mer (length of dna substring)
    Return data:
    ------------
        List of motifs which are presented in all dna strings with at most 2d hamming distance
    """
    patterns = []
    for i in range(0, len(dna[0]) - k + 1):
        k_mer = dna[0][i : i + k]
        k_mer_list = neighbours(k_mer, d)
        for k_mer in k_mer_list:
            for dna_str in dna[1::]:
                flag = True
                for j in range(0, len(dna_str) - k + 1):
                    compare_k_mer = dna_str[j : j + k]
                    if hamming_str(k_mer, compare_k_mer) <= d:
                        break
                else:
                    flag = False
                if not flag:
                    break
            else:
                patterns.append(k_mer)
    return patterns


@njit()
def pat_string_dists(pat: str, strng_list: list[str]) -> int:
    pat_length = len(pat)
    distance: int = 0
    for strng in strng_list:
        curr_dist = np.iinfo(np.int64).max
        for i in range(0, len(strng) - pat_length + 1):
            curr_sbstrng = strng[i : i + pat_length]
            tmp_dist = hamming_str(pat, curr_sbstrng)
            if tmp_dist < curr_dist:
                curr_dist = tmp_dist
        distance += curr_dist
    return distance


@njit(parallel=True)
def compute_profile_mat(
    dna_list: list[str], pseudocounts=True
) -> npt.NDArray[np.float64]:
    """
    Compute profile matrix from the list of dna stringd
    """
    if pseudocounts:
        profile_mat = np.ones((4, len(dna_list[0])), dtype=np.float64)
    else:
        profile_mat = np.zeros((4, len(dna_list[0])), dtype=np.float64)
    num_of_str = len(dna_list)
    dna_dct = {"A": 0, "C": 1, "G": 2, "T": 3}
    for dna in dna_list:
        # for idx, ncltd in enumerate(dna):
        for idx in prange(0, len(dna)):
            ncltd = dna[idx]
            profile_mat[dna_dct[ncltd], idx] += 1 / num_of_str
    return profile_mat


@njit(parallel=True)
def compute_motif_count_mat(dna_list: list[str]) -> npt.NDArray[np.int64]:
    """
    Compute count motifs matrix from the list of dna strings
    """
    profile_mat = np.zeros((4, len(dna_list[0])), dtype=np.int64)
    dna_dct = {"A": 0, "C": 1, "G": 2, "T": 3}
    for dna in dna_list:
        # for idx, ncltd in enumerate(dna):
        for idx in prange(0, len(dna)):
            ncltd = dna[idx]
            profile_mat[dna_dct[ncltd], idx] += 1
    return profile_mat


# @njit()
def score_motif_count_mat(motif_mat: npt.NDArray[np.int64]) -> int:
    """
    Compute score for count motifs matrix
    """
    return np.sum(motif_mat) - np.sum(np.max(motif_mat, axis=0))


@njit(parallel=True)
def compute_entropy(profile_mat: npt.NDArray[np.float64], eps: float = 1e-22) -> float:
    """
    Compute entropy of profile matrix
    Entropy of profile matrix is defined as sum of entropy of columns
    """
    return -np.sum(profile_mat * np.log2(profile_mat + eps))


def median_string_dna(dna_list: list[str], k: int) -> str:
    init_str = "A" * k
    return __median_string_dna__(dna_list, k, init_str)


@njit()
def __median_string_dna__(dna_list: list[str], k: int, init_str: str) -> str:
    """
    Find a motif which minimizes hamming distance score for motif matrix
    Inpit data:
    -----------
        dna_list - list of dna strings
        k - size of k-mers
        init_string - initial k-mer pattern, Needed for numba
    Output data:
    -----------
        Found motiff
    """
    fin_dist = np.iinfo(np.int64).max
    fin_pat = ""
    for pat in neighbours(init_str, k):
        curr_dist = pat_string_dists(pat, dna_list)
        if curr_dist < fin_dist:
            fin_dist = curr_dist
            fin_pat = pat
    return fin_pat


@njit()
def most_probable_k_mer(profile_mat: npt.NDArray[np.float64], dna: str, k: int) -> str:
    """
    Find most probable k-mer in dna string.
    Input data:
    -----------
        profile_mat - matrix with probabilities of each nucleotide in each position of dna string
        dna - dna string
        k - k-mer string size
    Output data:
    ------------
        Most probable k-mer which is presented in dna
    """
    dna_dct = {"A": 0, "C": 1, "G": 2, "T": 3}
    prb = 0
    fin_k_mer = dna[0 : 0 + k]
    for i in range(0, len(dna) - k + 1):
        curr_prb = 1
        k_mer = dna[i : i + k]
        for idx, ncltd in enumerate(k_mer):
            curr_prb *= profile_mat[dna_dct[ncltd], idx]
        if curr_prb > prb:
            prb = curr_prb
            fin_k_mer = k_mer
    return fin_k_mer


# @njit()
def greedy_motif_search(dna: list[str], k: int, t: int, pseudocounts=True) -> list[str]:
    """
    Greedy algorithm to search fot motif in list of dna strings
    Input data
    ----------
        dna - list of different dna strings or substrings
        k - size of k-mer
    Output data:
        Best found motif
    """
    dna_str_quant = t  # len(dna)
    best_motifs = [dna_str[0 : 0 + k] for dna_str in dna]
    # best_score = compute_entropy(compute_profile_mat(best_motifs))
    best_score = score_motif_count_mat(compute_motif_count_mat(List(best_motifs)))
    for i in range(0, len(dna[0]) - k + 1):
        curr_motif = [dna[0][i : i + k]]
        for j in range(1, dna_str_quant):
            curr_motif.append(
                most_probable_k_mer(
                    compute_profile_mat(List(curr_motif), pseudocounts=pseudocounts),
                    dna[j],
                    k,
                )
            )
        curr_score = score_motif_count_mat(compute_motif_count_mat(List(curr_motif)))
        if curr_score < best_score:
            best_motifs = curr_motif
            best_score = curr_score
    return best_motifs


@njit()
def randomized_motif_search(dna_list: list[str], k: int) -> tuple[list[str], float]:
    """
    Greedy, randomized algorithm to search fot motif in list of dna strings
    Input data
    ----------
        dna - list of different dna strings or substrings
        k - size of k-mer
    Output data:
        Best found motif
    """
    dna_quant = len(dna_list)
    start_idx = np.random.randint(low=0, high=len(dna_list[0]) - k + 1, size=dna_quant)
    motif = [dna_str[i : i + k] for dna_str, i in zip(dna_list, start_idx)]
    # score = score_motif_count_mat(compute_motif_count_mat(List(motif)))
    score = compute_entropy(compute_profile_mat(List(motif), pseudocounts=False))
    while True:
        curr_motif = []
        for j in range(0, dna_quant):
            curr_motif.append(
                most_probable_k_mer(
                    compute_profile_mat(List(motif), pseudocounts=True),
                    dna_list[j],
                    k,
                )
            )
        # curr_score = score_motif_count_mat(compute_motif_count_mat(List(curr_motif)))
        curr_score = compute_entropy(
            compute_profile_mat(List(curr_motif), pseudocounts=False)
        )
        if curr_score < score:
            score = curr_score
            motif = curr_motif
        else:
            return (motif, score)


@njit(parallel=True)
def gibbs_sampling_motif_finder(
    dna_list: list[str], k: int, N
) -> tuple[list[str], float]:
    """
    Greedy, randomized algorithm to search fot motif in list of dna strings
    Input data
    ----------
        dna - list of different dna strings or substrings
        k - size of k-mer
    Output data:
        Best found motif
    """
    dna_quant = len(dna_list)
    start_idx = np.random.randint(low=0, high=len(dna_list[0]) - k + 1, size=dna_quant)
    motif = [dna_str[i : i + k] for dna_str, i in zip(dna_list, start_idx)]
    best_motif = motif.copy()
    # score = score_motif_count_mat(compute_motif_count_mat(List(motif)))
    score = compute_entropy(compute_profile_mat(List(motif), pseudocounts=False))
    for _ in range(N):
        change_idx = np.random.randint(low=0, high=dna_quant)
        use_motif = List(motif[0:change_idx] + motif[change_idx + 1 : :])
        motif[change_idx] = most_probable_k_mer(
            compute_profile_mat(use_motif, pseudocounts=True),
            dna_list[change_idx],
            k,
        )
        # curr_score = score_motif_count_mat(compute_motif_count_mat(List(curr_motif)))
        curr_score = compute_entropy(compute_profile_mat(List(motif), pseudocounts=False))
        if curr_score < score:
            score = curr_score
            best_motif = motif.copy()
    return (best_motif, score)


if __name__ == "__main__":
    with open("input.txt", "r") as f:
        pat_str = f.readline().strip()
        dna_list = List(f.readline().strip().split())
    res = pat_string_dists(pat_str, dna_list)

    with open("output.txt", "w") as f:
        f.write(str(res))

    # with open("input.txt", "r") as f:
    #    k, d = [int(x) for x in f.readline().strip().split()]
    #    dna_list = List(f.readline().strip().split())
    # res = set(motif_enum(dna_list, k, d))

    # with open("output.txt", "w") as f:
    #    f.write(" ".join(res))
