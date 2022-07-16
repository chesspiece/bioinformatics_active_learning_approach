from three.grph_lib import deBrujin_grph

if __name__ == "__main__":
    with open("input_3_5.txt", "r") as f:
        dna_list = f.readline().strip().split()

    res = deBrujin_grph(dna_list)
    res.write_edges("output.txt")
