from typing import Optional

from consts.consts import CODON_TABLE, COMPLIMENTS, INTEGER_MASS
from one.reverse_complement import reversed_complement


def get_complement(strand: str, reverse: bool = True) -> str:
    """
    Return the reverse complement of the DNA sequence.

    :param strand: A sequence of nucleotides
    :param reverse: True if the complement should be reversed
    :return: The (reversed) complement of a DNA sequence
    """
    complement = "".join(COMPLIMENTS[x] for x in strand)
    if reverse:
        complement = complement[::-1]

    return complement


def translate_rna(
    rna: str, codon_table: dict[str, str] = CODON_TABLE, step_size: int = 3
) -> str:
    """
    Translate rna sequence into protein sequence
    Input data:
    -----------
        rna - input rna sequence
        cordon_table
        step_size
    Output data:
    ------------
        Aminoacid sequence
    """
    rna = rna.upper()
    aminoacid = []
    for i in range(0, len(rna) - step_size + 1, step_size):
        curr_cordon = rna[i : i + step_size]
        aminoacid.append(codon_table[curr_cordon])
    return "".join(aminoacid)


def dna2rna(dna: str) -> str:
    """
    Transcripe dna to rna
    """
    dna = dna.upper()
    return dna.replace("T", "U")


def find_peptide_encoding(dna_strand: str, peptide: str) -> list:
    """
    Given a peptide and a strand, finds the substrings that transcript in to RNA an then translates into the peptide.
    Since a DNA strand is a double-helix, the reverse compliment of each substring is also checked.

    :param strand: A sequence of nucleotides
    :param peptide: A sequence of amino acids
    :return: A list of substrings that translate into the peptide
    """
    str_len = 3 * len(peptide)
    substrings = []
    for i in range(len(dna_strand) - str_len + 1):
        substring = dna_strand[i : i + str_len]
        complement = reversed_complement(substring)

        rna = dna2rna(substring)
        if translate_rna(rna) == peptide:
            substrings.append(substring)
            continue
        rna = dna2rna(complement)
        if translate_rna(rna) == peptide:
            substrings.append(substring)
    return substrings


def linear_spectrum(
    peptide: str,
    amino_acid_mass: dict[str, int] = INTEGER_MASS,
    alphabet: Optional[list[str]] = None,
) -> list[int]:
    if alphabet is None:
        alphabet = list(peptide)
    prefix_mass = [0]
    for i in range(0, len(peptide)):
        prefix_mass.append(prefix_mass[i] + amino_acid_mass[peptide[i]])
    linear_spectrum = [0]
    for i in range(0, len(peptide)):
        for j in range(i + 1, len(peptide) + 1):
            linear_spectrum.append(prefix_mass[j] - prefix_mass[i])
    return sorted(linear_spectrum)


def cyclic_spectrum(
    peptide: str, amino_acid_mass=INTEGER_MASS, alphabet: Optional[list[str]] = None
) -> list[int]:
    if alphabet is None:
        alphabet = list(peptide)
    prefix_mass = [0]
    for i in range(0, len(peptide)):
        prefix_mass.append(prefix_mass[i] + amino_acid_mass[peptide[i]])
    peptide_mass = prefix_mass[-1]
    cyclic_spectrum = [0]
    for i in range(0, len(peptide)):
        for j in range(i + 1, len(peptide) + 1):
            cyclic_spectrum.append(prefix_mass[j] - prefix_mass[i])
            if (i > 0) and (j < (len(peptide))):
                cyclic_spectrum.append(peptide_mass - (prefix_mass[j] - prefix_mass[i]))
    return sorted(cyclic_spectrum)
