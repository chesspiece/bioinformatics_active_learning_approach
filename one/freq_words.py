# from collections import defaultdict
from numba import njit, prange


@njit()
def frequent_words(text: str, k: int) -> dict[str, int]:
    """
    Return dictionary with frequency of all k-mers in text
    Input data
    -----------
        k - length of pattern
        text - text pattern (usually dna)
    Output data
    -----------
        Dictionary with frequency of all k-mers in text
    """
    patt_dct: dict[str, int] = dict()
    for i in range(0, len(text) - k + 1):
        # for i in prange(0, len(text) - k):
        if text[i: i + k] in patt_dct:
            patt_dct[text[i: i + k]] += 1
        else:
            patt_dct[text[i: i + k]] = 1
    return patt_dct


if __name__ == "__main__":
    with open("input.txt", "r") as f:
        dna_str = f.readline().strip()
        k = int(f.readline().strip())
    ret = frequent_words(dna_str, k)

    most_frequent_kmers: list[str] = []
    v_prev = 0
    for k, v in sorted(ret.items(), key=lambda x: x[1], reverse=True):
        if v_prev > v:
            break
        most_frequent_kmers.append(k)
        v_prev = v

    with open("output.txt", "w") as f:
        f.write(" ".join(most_frequent_kmers))
