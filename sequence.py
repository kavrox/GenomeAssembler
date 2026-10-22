"""Part 1 - String ADT for DNA sequences.

A DNA sequence is modelled as an immutable Python ``str`` over the alphabet
{A, C, G, T}. Python strings are array-backed, so ``len`` and indexing are
O(1) and slicing a substring of length m costs O(m) (it copies).
"""

from __future__ import annotations

DNA_ALPHABET = frozenset("ACGT")
_COMPLEMENT_TABLE = str.maketrans("ACGT", "TGCA")


def validate_dna(seq: str, name: str | None = None) -> str:
    """Normalise a DNA string (drop whitespace, uppercase) and check the alphabet.

    Raises ValueError for an empty sequence or any character outside A/C/G/T
    (including ambiguity codes such as N). O(n).
    """
    if seq is None:
        raise ValueError("sequence is missing")
    cleaned = "".join(str(seq).split()).upper()
    where = f" in {name!r}" if name else ""
    if not cleaned:
        raise ValueError(f"empty sequence{where}")
    bad = set(cleaned) - DNA_ALPHABET
    if bad:
        raise ValueError(
            f"invalid DNA character(s) {sorted(bad)}{where}; only A, C, G, T are allowed"
        )
    return cleaned


def is_valid_dna(seq: str) -> bool:
    try:
        validate_dna(seq)
        return True
    except ValueError:
        return False


def complement(seq: str) -> str:
    """A<->T, C<->G for every base. O(n)."""
    return seq.translate(_COMPLEMENT_TABLE)


def reverse(seq: str) -> str:
    return seq[::-1]


def reverse_complement(seq: str) -> str:
    """Sequence of the opposite strand, read 5'->3'. ATGC -> GCAT. O(n)."""
    return seq.translate(_COMPLEMENT_TABLE)[::-1]


def canonical(kmer: str) -> str:
    """Lexicographically smaller of a k-mer and its reverse complement. O(k)."""
    rc = reverse_complement(kmer)
    return min(kmer, rc)


def gc_content(seq: str) -> float:
    """Fraction of G and C bases (0.0 for an empty string)."""
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def prefix(seq: str, m: int) -> str:
    return seq[:m]


def suffix(seq: str, m: int) -> str:
    return seq[len(seq) - m:] if m > 0 else ""
