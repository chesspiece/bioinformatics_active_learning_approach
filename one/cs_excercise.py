from numba import njit
from .skew import hamming_str


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
    suffix_neighboors = neighbours(dna_pat[1::], d=d)
    for subpat in suffix_neighboors:
        if hamming_str(dna_pat[1::], subpat) < d:
            # if hs_jit(dna_pat[1::], subpat) < d:
            for ncl_idx in ["A", "G", "C", "T"]:
                neighborhood.append(ncl_idx + subpat)
        else:
            neighborhood.append(dna_pat[0] + subpat)
    return neighborhood


if __name__ == "__main__":
    with open("input.txt", "r") as f:
        dna_motif = f.readline().strip()
        d = int(f.readline().strip())

    nghbrhd = neighbours(dna_motif, d)
    with open("output.txt", "w") as f:
        f.write(" ".join(nghbrhd))
