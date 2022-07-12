import numpy as np
import numpy.typing as npt
from numba import njit
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


@njit()
def compute_profile_mat(dna_list: list[str]) -> npt.NDArray[np.float64]:
    """
    Compute profile matrix from the list of dna stringd
    """
    profile_mat = np.zeros((4, len(dna_list[0])), dtype=np.float64)
    num_of_str = len(dna_list)
    dna_dct = {"A": 0, "C": 1, "G": 2, "T": 3}
    for dna in dna_list:
        for idx, ncltd in enumerate(dna):
            profile_mat[dna_dct[ncltd], idx] += 1 / num_of_str
    return profile_mat


def compute_entropy(profile_mat: npt.NDArray[np.float64]) -> float:
    """
    Compute entropy of profile matrix
    Entropy of profile matrix is defined as sum of entropy of columns
    """
    # tmp = profile_mat * np.log2(profile_mat)
    tmp = profile_mat * np.where(np.isclose(profile_mat, 0), 0, np.log2(profile_mat))
    # profile_mat * np.where(np.isclose(profile_mat, 0), 0, np.log2(profile_mat))
    # tmp[np.isnan(tmp)] = 0
    return -np.sum(tmp)


@njit()
def median_string_dna(dna_list: list[str], k: int, init_str: str) -> str:
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
