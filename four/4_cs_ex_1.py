from four.four_lib import cyclic_spectrum

if __name__ == "__main__":
    with open("input_4_cs_1.txt") as f:
        peptide = f.readline().strip()
    answ = cyclic_spectrum(peptide)

    with open("output.txt", "w") as f:
        f.write(" ".join(map(str, answ)))
