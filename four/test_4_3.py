from four.four_lib import leaderboard_sequencing, mass

if __name__ == "__main__":
    with open("input_4_5.txt") as f:
        N = int(f.readline().strip())
        spectrum = [int(x) for x in f.readline().strip().split()]
    leader_peptide = leaderboard_sequencing(spectrum, N)

    answ = [str(mass(i)) for i in leader_peptide]
    with open("output.txt", "w") as f:
        f.write("-".join(map(str, answ)))
