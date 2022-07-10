import numpy as np
import numpy.typing as npt

from numba import njit


def convert_DNA(dna: str, conv: dict[str, int]) -> npt.NDArray[np.int8]:
    dna_number_arr = np.zeros((len(dna)), dtype=np.int8)
    for idx, nucleotide in enumerate(dna):
        dna_number_arr[idx] = conv[nucleotide]
    return dna_number_arr


def compute_skew(
    dna_arr: npt.NDArray[np.int8], convert: dict[str, int]
) -> npt.NDArray[np.int64]:
    g_n = convert["G"]
    c_n = convert["C"]
    skew = np.zeros(dna_arr.size + 1, dtype=np.int64)
    skew[1:] = dna_arr.copy()
    skew[(skew != g_n) & (skew != c_n)] = 0
    skew[skew == g_n] = 1
    skew[skew == c_n] = -1
    return skew.cumsum()


@njit()
def hamming_str(dna: str, dna_compare: str) -> int:
    return sum([0 if ncl_1 == ncl_2 else 1 for ncl_1, ncl_2 in zip(dna, dna_compare)])
    #return np.where(
    #    [ncl_1 == ncl_2 for ncl_1, ncl_2 in zip(dna, dna_compare)], 0, 1
    #).sum()


def hamming(
    dna_arr: npt.NDArray[np.int8], dna_arr_compare: npt.NDArray[np.int8]
) -> int:
    return np.where(dna_arr == dna_arr_compare, 0, 1).sum()


def hamming_idx(
    dna_arr: npt.NDArray[np.int8], dna_pattern: npt.NDArray[np.int8], k: int
) -> list[int]:
    res: list[int] = []
    for i in range(dna_arr.size - dna_pattern.size + 1):
        if hamming(dna_arr[i : i + dna_pattern.size], dna_pattern) <= k:
            res.append(i)
    return res


def count_pat_mismatch(
    dna_arr: npt.NDArray[np.int8], dna_pattern: npt.NDArray[np.int8], k: int
) -> int:
    cnt = 0
    for i in range(dna_arr.size - dna_pattern.size + 1):
        if hamming(dna_arr[i : i + dna_pattern.size], dna_pattern) <= k:
            cnt += 1
    return cnt


if __name__ == "__main__":
    """
    with open("input.txt", "r") as f:
        dna_str = f.readline().strip()
        k, L, t = map(int, f.readline().strip().split())
    answ = set(find_clumps(dna_str, k, L, t))
    with open("output.txt", "w") as f:
        f.write(" ".join(answ))
    """

    with open("one/input.txt", "r") as f:
        dna_str_comp = f.readline().strip()
        dna_str = f.readline().strip()
        k = int(f.readline().strip())
    conv = {"A": 1, "C": 2, "G": 3, "T": 4, "N": 0}
    conv_rev = {v: k for k, v in conv.items()}
    dna_arr = convert_DNA(dna_str, conv)
    dna_arr_comp = convert_DNA(dna_str_comp, conv)
    # skew = compute_skew(dna_arr, conv)
    # min_skew_idx = np.where(skew == skew.min())[0]

    # hd = hamming_idx(dna_arr, dna_arr_comp, k)
    hd = count_pat_mismatch(dna_arr, dna_arr_comp, k)

    with open("output.txt", "w") as f:
        # f.write(" ".join(map(str, hd)))
        f.write(str(hd))
