if __name__ == "__main__":
    with open("input_3_1.txt", "r") as f:
        k = int(f.readline().strip())
        dna = f.readline().strip()
    res = []
    for i in range(0, len(dna) - k + 1):
        res.append(dna[i : i + k])

    with open("output.txt", "w") as f:
        f.write(" ".join(res))
