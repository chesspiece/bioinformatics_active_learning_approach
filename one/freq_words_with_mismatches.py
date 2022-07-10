from cs_excercise import neighbours
from numba import njit
from reverse_complement import reversed_complement as fnc


@njit()
def frequent_words_w_m(text: str, k: int, d: int) -> dict[str, int]:
    """
    Return dictionary with frequency of all k-mers in text
    Input data:
    -----------
        d - maximum alowable msimatches
        k - length of pattern
        text - text pattern (usually dna)
    Output data:
    -----------
        Dictionary with frequency of all k-mers in text
    """
    # fnc = njit(reversed_complement)
    patt_dct: dict[str, int] = dict()
    for i in range(0, len(text) - k):
        pat = text[i : i + k]
        pat_neigh = neighbours(pat, d) + neighbours(fnc(pat), d)
        for i in range(len(pat_neigh)):
            nbr_pat = pat_neigh[i]
            if nbr_pat in patt_dct:
                patt_dct[nbr_pat] += 1
            else:
                patt_dct[nbr_pat] = 1
    return patt_dct


if __name__ == "__main__":
    with open("input.txt", "r") as f:
        dna_motif = f.readline().strip()
        k, d = map(int, f.readline().strip().split())

    ret = frequent_words_w_m(dna_motif, k, d)

    most_frequent_kmers: list[str] = []
    v_prev = 0
    for k, v in sorted(ret.items(), key=lambda x: x[1], reverse=True):
        if v_prev > v:
            break
        most_frequent_kmers.append(k)
        v_prev = v

    with open("output.txt", "w") as f:
        f.write(" ".join(most_frequent_kmers))
