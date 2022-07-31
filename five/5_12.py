from five.dp_lib import middle_node

if __name__ == "__main__":
    with open("input_5_12.txt") as f:
        match_rewad, mismatch_penalty, indel_penalty = [
            int(x) for x in f.readline().strip().split()
        ]
        dna_1 = f.readline().strip()
        dna_2 = f.readline().strip()

    (i, j, s_i, s_j), _ = middle_node(dna_2, dna_1, match_rewad, mismatch_penalty, indel_penalty)

    with open("output.txt", "w") as f:
       f.write(" ".join(map(str, [i, j])) + "\n") 
       f.write(" ".join(map(str, [s_i, s_j])) + "\n") 
