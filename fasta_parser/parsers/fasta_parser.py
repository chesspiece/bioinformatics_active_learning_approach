# -*- coding: utf-8 -*-
"""
Implementation of FASTA format parser
"""
from __future__ import annotations

from typing import Iterator, List, TextIO, Tuple


def yield_line(f_d: TextIO) -> Iterator[Tuple[str, str]]:
    """
    Generator function to read Fasta files
    Input:
    ------
        f_d - descriptor of text file
    Output:
    -------
        Tuple of FASTA descriptor and full DNA sequence
    """
    for line in f_d:
        if line[0] == ">":
            title: str = line[1:].strip()
            break
    else:
        return
    dna_sbln: List[str] = []
    for line in f_d:
        if line[0] == ">":
            yield title, "".join(dna_sbln).replace(" ", "").replace("\r", "")
            title = line[1:].strip()
            dna_sbln = []
            continue
        dna_sbln.append(line.strip())
    yield title, "".join(dna_sbln).replace(" ", "").replace("\r", "")


class FastaParser:
    """
    Calss to pasrse FASTA format
    """

    def __init__(self, fname: str, enc: str = "utf-8") -> None:
        """
        Initialize file name and file mode
        """
        self.fname: str = fname
        # self.mode: str = mode
        self.f_descr: TextIO
        self.seq_it: Iterator[Tuple[str, str]]
        self.encoding: str = enc

    def __iter__(self) -> Iterator[Tuple[str, str]]:
        self.f_descr.seek(0)
        return yield_line(self.f_descr)

    def open(self) -> None:
        """
        Open FASTA file connected to this class
        """
        # should add tests!
        # pylint: disable-next=consider-using-with
        self.f_descr = open(self.fname, "r", encoding=self.encoding)
        # self.seq_it = yield_line(self.f_descr)

    def close(self) -> None:
        """
        Close FASTA file connected to this class
        """
        self.f_descr.close()

    def __enter__(self) -> FastaParser:
        """
        Custom context manager "enter"
        """
        self.open()
        return self

    # pylint: disable-next=redefined-builtin
    def __exit__(self, type, value, traceback) -> None:
        """
        Custom context manager "exit"
        """
        self.close()
