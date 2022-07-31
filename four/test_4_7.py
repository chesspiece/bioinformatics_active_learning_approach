from four.four_lib import leaderboard_sequencing_spectrum, mass
from consts.consts import EXTENDED_AMINO

if __name__ == "__main__":
    with open("input_4_7.txt") as f:
        M = int(f.readline().strip())
        N = int(f.readline().strip())
        spectrum = [int(x) for x in f.readline().strip().split()]
    leader_peptide, alphabet = leaderboard_sequencing_spectrum(spectrum, M, N)

    answ = [str(mass(i, mass_table=EXTENDED_AMINO)) for i in leader_peptide]
    with open("output.txt", "w") as f:
        f.write("-".join(map(str, answ)))
