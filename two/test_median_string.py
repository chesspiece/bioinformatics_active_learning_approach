from numba.typed import List

from two.motif_lib import median_string_dna

if __name__ == "__main__":
    with open("input.txt", "r") as f:
        k = int(f.readline().strip())
        dna_list = List(f.readline().strip().split())
    res = median_string_dna(dna_list, k, "A" * k)

    with open("output.txt", "w") as f:
        f.write(res)
