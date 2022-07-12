from two.motif_lib import compute_profile_mat, compute_entropy
from numba.typed import List

if __name__ == "__main__":
    dna_list: list[str] = []
    with open("input.txt", "r") as f:
        for line in f:
            dna_list.append(line.strip())
    dna_list = List(dna_list)
    res = compute_entropy(compute_profile_mat(dna_list))
    print(res)
