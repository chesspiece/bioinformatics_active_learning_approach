from five.dp_lib import output_lcs

if __name__ == "__main__":
    with open("input_5_8.txt") as f:
        dna_1 = f.readline().strip()
        dna_2 = f.readline().strip()

    score, outp1, outp2 = output_lcs(
        dna_1, dna_2
    )

    with open("output.txt", "w") as f:
        f.write(str(score) + "\n")
        f.write("\n".join([outp1, outp2]))
