"""Part 1 - k-mer generation (sliding window) and counting (hash map).

A string of length L has exactly L - k + 1 k-mers: the window can start at
positions 0, 1, ..., L - k. Counting uses a hash map (collections.Counter, a
dict subclass) so every insert/lookup is O(1) expected.
"""

from __future__ import annotations

from collections import Counter


def kmer_count_for_length(length: int, k: int) -> int:
    """Number of k-mers in a string of the given length: max(0, L - k + 1)."""
    return max(0, length - k + 1)


def iter_kmers(seq: str, k: int):
    """Yield k-mers left to right. Each slice copies k characters -> O(L*k) total."""
    if k < 1:
        raise ValueError("k must be >= 1")
    for i in range(len(seq) - k + 1):
        yield seq[i:i + k]


def generate_kmers(seq: str, k: int) -> list:
    """All L - k + 1 overlapping k-mers of ``seq`` ([] when k > L)."""
    return list(iter_kmers(seq, k))


def count_kmers(reads, k: int) -> Counter:
    """Hash map {k-mer: number of times observed across all reads} (= coverage).

    Time  O(K * k) where K = sum over reads of (L - k + 1)  (hashing a k-mer is O(k))
    Space O(D * k) where D = number of distinct k-mers
    """
    if k < 1:
        raise ValueError("k must be >= 1")
    freq = Counter()
    for read in reads:
        freq.update(iter_kmers(read, k))
    return freq


def kmer_spectrum(freq: Counter) -> dict:
    """Histogram {coverage: how many distinct k-mers have that coverage}.

    Error k-mers pile up at coverage 1-2, real k-mers near the sequencing depth,
    which is how a sensible ``min_coverage`` threshold is chosen.
    """
    return dict(sorted(Counter(freq.values()).items()))


def total_kmers(reads, k: int) -> int:
    return sum(kmer_count_for_length(len(r), k) for r in reads)
