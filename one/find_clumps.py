from freq_words import frequent_words
from freq_words_with_mismatches import frequent_words_w_m
from numba import njit

from fasta_parser.parsers.fasta_parser import FastaParser


@njit()
def find_clumps(text: str, k: int, L: int, t: int, d: int = 3):
    frequent_patterns: list[str] = []
    #freq_table = frequent_words(text[0 : 0 + L], k)
    freq_table = frequent_words_w_m(text[0 : 0 + L], k, d)
    for pat_name, freq in sorted(freq_table.items(), key=lambda x: x[1], reverse=True):
        if freq >= t:
            frequent_patterns.append(pat_name)

    for i in range(1, len(text) - L):
        str_ret = text[i - 1 : i - 1 + k]
        str_add = text[i + L - k - 1 : i + L - 1]
        freq_table[str_ret] -= 1
        if str_add in freq_table:
            freq_table[str_add] += 1
        else:
            freq_table[str_add] = 1
        pat_name, freq = (str_add, freq_table[str_add])
        if freq >= t:
            frequent_patterns.append(pat_name)
    return frequent_patterns


if __name__ == "__main__":
    """
    with open("input.txt", "r") as f:
        dna_str = f.readline().strip()
        k, L, t = map(int, f.readline().strip().split())
    answ = set(find_clumps(dna_str, k, L, t))
    with open("output.txt", "w") as f:
        f.write(" ".join(answ))
    """

    #with open("Salmonella_enterica.txt", "r") as f:
    #    dna_str = f.readline().strip()
    dna_str = ""
    with FastaParser("Salmonella_enterica.txt") as fh:
        for _, dna in fh:
            dna_str = dna
    dna_str = dna_str[3764856-1000:3764956+1000]

    k, L, t = (9, 500, 3)
    answ = set(find_clumps(dna_str, k, L, t))
    with open("output.txt", "w") as f:
        f.write(" ".join(answ))
