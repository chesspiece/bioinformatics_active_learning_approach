from three.grph_lib import deBrujin_grph, str_comp

if __name__ == "__main__":
    with open("input_3_4.txt", "r") as f:
        k = int(f.readline().strip())
        dna_list = str_comp(f.readline().strip(), k)

    res = deBrujin_grph(dna_list)
    res.write_edges("output.txt")
