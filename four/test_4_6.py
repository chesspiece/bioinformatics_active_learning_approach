from four.four_lib import spectral_convolution

if __name__ == "__main__":
    with open("input_4_6.txt") as f:
        spectrum = [int(x) for x in f.readline().strip().split()]
    convl = spectral_convolution(spectrum)

    with open("output.txt", "w") as f:
        f.write(" ".join(map(str, convl)))
