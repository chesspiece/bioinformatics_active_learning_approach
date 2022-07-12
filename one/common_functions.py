from numba import njit

from .skew import hamming_str


@njit()
def neighbours(
    dna_pat: str,
    d=0,
) -> list[str]:
    """
    Generate d-neighbourhood of dna_pat.
    d-neighborhood is a set of all strings
    which have at most d hamming distance
    from dna_pat

    Input data:
    -----------
        dna_pat - inpu pattern
        d - maximum allowable hamming distance in generated neighboorhood strings
    """
    if d == 0:
        return [dna_pat]
    if len(dna_pat) == 1:
        return ["A", "G", "C", "T"]
    neighborhood: list[str] = []
    suffix_neighboors = neighbours(dna_pat[1::], d=d)
    for subpat in suffix_neighboors:
        if hamming_str(dna_pat[1::], subpat) < d:
            for ncl_idx in ["A", "G", "C", "T"]:
                neighborhood.append(ncl_idx + subpat)
        else:
            neighborhood.append(dna_pat[0] + subpat)
    return neighborhood


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
    for i in range(0, len(text) - k):
        if text[i : i + k] in patt_dct:
            patt_dct[text[i : i + k]] += 1
        else:
            patt_dct[text[i : i + k]] = 1
    return patt_dct
