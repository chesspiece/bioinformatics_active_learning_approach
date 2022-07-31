from five.dp_lib import output_lcs

if __name__ == "__main__":
    with open("input_5_11.txt") as f:
        match_rewad, mismatch_penalty, indel_penalty = [
            int(x) for x in f.readline().strip().split()
        ]
        dna_1 = f.readline().strip()
        dna_2 = f.readline().strip()

    score, outp1, outp2 = output_lcs(
        dna_1, dna_2, match_rewad, mismatch_penalty, indel_penalty
    )

    with open("output.txt", "w") as f:
        f.write(str(score) + "\n")
        f.write("\n".join([outp1, outp2]))
