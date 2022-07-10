from one.common_functions import neighbours
from one.skew import hamming_str


def motif_enum(dna: list[str], k: int, d: int) -> set[str]:
    """
    Find if motif presented in all dna strings or substrings with at most 2d mismatches
    Input data:
    -----------
        dna - input dna string
        d - allowed number of mismathces
        k - k-mer (length of dna substring)
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
    return set(patterns)


if __name__ == "__main__":
    with open("input.txt", "r") as f:
        k, d = [int(x) for x in f.readline().strip().split()]
        dna_list = f.readline().strip().split()
    res = motif_enum(dna_list, k, d)

    with open("output.txt", "w") as f:
        f.write(" ".join(res))
