"""Part 3 - Synthetic data: random genomes, sequencing reads and mutations.

Used for bigger demos, for the benchmark and for testing the full pipeline
against a known truth.
"""

from __future__ import annotations

import random

BASES = "ACGT"


def random_genome(length: int, seed: int = None) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice(BASES) for _ in range(length))


def simulate_reads(genome: str, read_length: int = 50, coverage: float = 10.0,
                   error_rate: float = 0.0, seed: int = None) -> list:
    """Error-free tiling (step = read_length // 2) plus random reads up to ``coverage``.

    The tiling guarantees that consecutive reads overlap by >= read_length // 2,
    so every k-mer with k <= read_length // 2 + 1 is observed at least once.
    Substitution errors are injected with probability ``error_rate`` per base.
    """
    if read_length > len(genome):
        return [genome]
    rng = random.Random(seed)
    g, rl = len(genome), read_length
    starts = list(range(0, g - rl + 1, max(1, rl // 2)))
    if starts[-1] != g - rl:
        starts.append(g - rl)
    needed = int(coverage * g / rl)
    while len(starts) < needed:
        starts.append(rng.randint(0, g - rl))
    reads = []
    for s in starts:
        read = list(genome[s:s + rl])
        if error_rate > 0:
            for i, base in enumerate(read):
                if rng.random() < error_rate:
                    read[i] = rng.choice([b for b in BASES if b != base])
        reads.append("".join(read))
    rng.shuffle(reads)
    return reads


def _make_event(seq: str, p: int, kind: str, rng):
    """Return (ref allele, alt allele) for an event at 0-based index p, or None if the
    event would be ambiguous (an indel that could be shifted left/right in a run)."""
    if kind == "SUBSTITUTION":
        return seq[p], rng.choice([b for b in BASES if b != seq[p]])
    n = rng.randint(1, 3)
    if kind == "INSERTION":
        ins = "".join(rng.choice(BASES) for _ in range(n))
        if ins[0] == seq[p] or ins[-1] == seq[p - 1]:
            return None
        return "", ins
    if seq[p] == seq[p + n] or seq[p - 1] == seq[p + n - 1]:
        return None
    return seq[p:p + n], ""


def mutate(seq: str, n_snv: int = 1, n_ins: int = 0, n_del: int = 0,
           min_gap: int = 30, seed: int = None):
    """Apply well-separated SNVs and unambiguous 1-3 bp insertions / deletions.

    Returns (mutated sequence, list of (type, 1-based ref position, ref, alt)),
    with the same position convention as the mutation analyzer (an insertion is
    reported before the given reference position). Events are at least
    ``min_gap`` bases apart so each one is resolvable by the anchor algorithm.
    """
    rng = random.Random(seed)
    kinds = ["SUBSTITUTION"] * n_snv + ["INSERTION"] * n_ins + ["DELETION"] * n_del
    rng.shuffle(kinds)
    events, tries = [], 0
    for kind in kinds:
        while True:
            tries += 1
            if tries > 100_000:
                raise ValueError("sequence too short for that many separated mutations")
            p = rng.randint(min_gap, len(seq) - min_gap - 4)
            if any(abs(p - q) < min_gap for q, _, _, _ in events):
                continue
            alleles = _make_event(seq, p, kind, rng)
            if alleles is not None:
                events.append((p, kind) + alleles)
                break
    events.sort()
    out, truth, cursor = [], [], 0
    for p, kind, ref, alt in events:
        out.append(seq[cursor:p])
        out.append(alt)
        cursor = p + len(ref)
        truth.append((kind, p + 1, ref, alt))
    out.append(seq[cursor:])
    return "".join(out), truth
