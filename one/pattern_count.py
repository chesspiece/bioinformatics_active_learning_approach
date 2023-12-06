import re
from typing import Pattern


def pattern_count(text: str, pattern: Pattern[str]) -> list[str]:
    return pattern.findall(text)


def pattern_index(text: str, pattern: Pattern[str]) -> list[int]:
    return [x.start() for x in pattern.finditer(text)]


if __name__ == "__main__":
    # with open("Vibrio_cholerae.txt", "r") as f:
    #    #pattern = f.readline().strip()
    #    dna_str = f.readline().strip()

    with open("./input/input_1_1.txt", "r") as f:
        dna_str = f.readline().strip()
        pattern = f.readline().strip()

    # pattern = "CTTGATCAT"
    re_pattern = re.compile("(?=(" + pattern + "))")
    with open("output.txt", "w") as f:
        # f.write(" ".join((map(str, pattern_index(dna_str, re_pattern)))))
        f.write(str(len(pattern_count(dna_str, re_pattern))))
