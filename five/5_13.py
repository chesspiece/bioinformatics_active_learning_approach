import sys
from five.dp_lib import linear_space_aligment, find_path_back

if __name__ == "__main__":
    sys.setrecursionlimit(150000)
    with open("input_5_13.txt") as f:
        match_rewad, mismatch_penalty, indel_penalty = [
            int(x) for x in f.readline().strip().split()
        ]
        dna_1 = f.readline().strip()
        dna_2 = f.readline().strip()

    path = []
    linear_space_aligment(
        dna_1,
        dna_2,
        0,
        len(dna_1),
        0,
        len(dna_2),
        path,
        match_rewad,
        mismatch_penalty,
        indel_penalty,
    )
    print("Yay")
    score, al1, al2 = find_path_back(dna_1, dna_2, path, match_rewad, mismatch_penalty, indel_penalty)

    with open("output.txt", "w") as f:
        f.write("\n".join([str(score), al1, al2]))
