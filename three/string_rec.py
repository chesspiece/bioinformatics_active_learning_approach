if __name__ == "__main__":
    with open("./input/input_3_2.txt", "r") as f:
        dna_list = f.readline().strip().split()
    res = dna_list[0]
    for i in range(1, len(dna_list)):
        res = res + dna_list[i][-1]

    with open("output.txt", "w") as f:
        f.write("".join(res))
