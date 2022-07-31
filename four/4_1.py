from four.four_lib import translate_rna

if __name__ == "__main__":
    with open("input_4_1.txt") as f:
        rna = f.readline().strip()
    answ = translate_rna(rna)
    with open("output.txt", "w") as f:
        f.write("".join(answ))
