from copy import deepcopy

from .skew import hamming_str
import numpy as np
import numpy.typing as npt
from numba import njit, prange

@njit()
def neighbours(
    dna_pat: str,
    d=0,
) -> list[str]:
    if d == 0:
        return [dna_pat]
    if len(dna_pat) == 1:
        return ["A", "G", "C", "T"]
    neighborhood: list[str] = []
    suffix_neighboors = neighbours(deepcopy(dna_pat[1::]), d=d)
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
        # for i in prange(0, len(text) - k):
        if text[i : i + k] in patt_dct:
            patt_dct[text[i : i + k]] += 1
        else:
            patt_dct[text[i : i + k]] = 1
    return patt_dct
