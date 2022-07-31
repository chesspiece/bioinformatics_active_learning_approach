from four.four_lib import score_peptide, cyclic_spectrum

if __name__ == "__main__":
    with open("input_4_4.txt") as f:
        peptide = f.readline().strip()
        full_spectrum = [int(x) for x in f.readline().strip().split()]
    test_spectrum = cyclic_spectrum(peptide)

    answ = score_peptide(test_spectrum, full_spectrum)

    with open("output.txt", "w") as f:
        f.write(str(answ))
