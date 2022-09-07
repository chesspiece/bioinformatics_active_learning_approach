from four.four_lib import find_peptide_encoding

if __name__ == "__main__":
    with open("input_4_2.txt") as f:
        dna = f.readline().strip()
        peptide_pat = f.readline().strip()
    # answ = find_peptide_encoding(dna, peptide_pat)
    answ = find_peptide_encoding(dna, peptide_pat)
    with open("output.txt", "w") as f:
        f.write("\n".join(answ))
