from three.grph_lib import myGrph

if __name__ == "__main__":
    with open("./input/input_3_3.txt", "r") as f:
        dna_list = f.readline().strip().split()
    res = myGrph(dna_list)
    res.write_edges("output.txt")
