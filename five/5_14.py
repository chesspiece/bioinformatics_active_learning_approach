from five.dp_lib import output_lcs_3d

if __name__ == "__main__":
    with open("input_5_14.txt") as f:
        dna_1 = f.readline().strip()
        dna_2 = f.readline().strip()
        dna_3 = f.readline().strip()

    path = []
    res = output_lcs_3d(dna_1, dna_2, dna_3)

    with open("output.txt", "w") as f:
        f.write("\n".join(map(str, res)))
