import pytest

from parsers.fasta_parser import FastaParser


def test_fasta_parser():
    expected_rslt = ["CTGACCCCGATCAAGG", "GATTTCATGTT", "CAATCGTCA"]
    archieved_rslt = []
    with FastaParser("./tests/test_fasta.txt") as fh:
        for _, dna in fh:
            archieved_rslt.append(dna)
    assert archieved_rslt == expected_rslt


if __name__ == "__main__":
    test_fasta_parser()
