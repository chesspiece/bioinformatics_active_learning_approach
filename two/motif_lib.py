import numpy as np
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
