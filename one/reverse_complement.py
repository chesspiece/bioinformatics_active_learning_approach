from numba import njit


@njit()
def reversed_complement(text: str) -> str:

    complement_dcit = {"A": "T", "T": "A", "G": "C", "C": "G"}
    answ: list[str] = []
    # for c in reversed(text):
    text_len = len(text)
    for i in range(text_len):
        c = text[text_len - i - 1]
        answ.append(complement_dcit[c])
    return "".join(answ)


if __name__ == "__main__":
    with open("input.txt", "r") as f:
        dna_str = f.readline().strip()
    with open("output.txt", "w") as f:
        f.write(reversed_complement(dna_str))
