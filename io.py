"""Part 1 - Input manager: plain text, FASTA and FASTQ parsing + FASTA output.

Every source handed to ``load_reads`` / ``load_sequence`` can be either a file
path or a literal DNA string, so the same command works for files and for test
cases typed during a live demo.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field

from .sequence import validate_dna


@dataclass
class SequenceRecord:
    name: str
    sequence: str
    quality: list = field(default=None, repr=False)  # Phred scores (FASTQ only)

    @property
    def mean_quality(self):
        if not self.quality:
            return None
        return sum(self.quality) / len(self.quality)


def detect_format(text: str) -> str:
    """'fasta' if the first non-blank line starts with '>', 'fastq' for '@', else 'plain'."""
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            return "fasta"
        if line.startswith("@"):
            return "fastq"
        return "plain"
    return "plain"


def parse_fasta(text: str) -> list:
    """Return SequenceRecords; multi-line sequences are joined and validated."""
    records, name, buf = [], None, []

    def flush():
        if name is not None:
            records.append(SequenceRecord(name, validate_dna("".join(buf), name)))

    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            flush()
            name, buf = (line[1:].strip() or f"seq{len(records) + 1}"), []
        else:
            if name is None:
                raise ValueError("FASTA sequence data found before the first '>' header")
            buf.append(line)
    flush()
    return records


def parse_fastq(text: str) -> list:
    """Return SequenceRecords with Phred+33 qualities (Q = ord(char) - 33)."""
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    if len(lines) % 4 != 0:
        raise ValueError(f"FASTQ must have 4 lines per read, got {len(lines)} non-blank lines")
    records = []
    for i in range(0, len(lines), 4):
        header, seq, plus, qual = lines[i:i + 4]
        lineno = i + 1
        if not header.startswith("@"):
            raise ValueError(f"FASTQ record at line {lineno}: header must start with '@'")
        if not plus.startswith("+"):
            raise ValueError(f"FASTQ record at line {lineno}: third line must start with '+'")
        name = header[1:].strip() or f"read{len(records) + 1}"
        if len(seq) != len(qual):
            raise ValueError(f"FASTQ read {name!r}: sequence and quality lengths differ")
        scores = [ord(c) - 33 for c in qual]
        if any(q < 0 or q > 93 for q in scores):
            raise ValueError(f"FASTQ read {name!r}: quality characters outside Phred+33 range")
        records.append(SequenceRecord(name, validate_dna(seq, name), scores))
    return records


def parse_plain(text: str) -> list:
    """One sequence per non-blank line; lines starting with '#' are comments."""
    records = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        name = f"seq{len(records) + 1}"
        records.append(SequenceRecord(name, validate_dna(line, name)))
    return records


def parse_text(text: str) -> list:
    fmt = detect_format(text)
    if fmt == "fasta":
        return parse_fasta(text)
    if fmt == "fastq":
        return parse_fastq(text)
    return parse_plain(text)


def filter_by_quality(records: list, min_quality: float) -> tuple:
    """Split records into (kept, dropped) by mean Phred quality. Records without
    qualities (FASTA / plain) are always kept."""
    kept, dropped = [], []
    for r in records:
        mq = r.mean_quality
        (dropped if mq is not None and mq < min_quality else kept).append(r)
    return kept, dropped


def read_records(path: str) -> list:
    with open(path, "r", encoding="utf-8") as fh:
        text = fh.read()
    try:
        records = parse_text(text)
    except ValueError as exc:
        raise ValueError(f"{path}: {exc}") from None
    if not records:
        raise ValueError(f"{path}: no sequences found")
    return records


def _looks_like_path(src: str) -> bool:
    return any(ch in src for ch in "./\\:")


def _literal(src: str, label: str) -> str:
    if _looks_like_path(src):
        raise ValueError(f"file not found: {src}")
    return validate_dna(src, label)


def load_reads(sources, min_quality: float = 0) -> tuple:
    """Load reads from file paths and/or literal sequences.

    Returns (records, dropped) where ``dropped`` are FASTQ reads whose mean
    quality is below ``min_quality``.
    """
    if isinstance(sources, str):
        sources = [sources]
    records = []
    for src in sources:
        if os.path.isfile(src):
            records.extend(read_records(src))
        else:
            name = f"read{len(records) + 1}"
            records.append(SequenceRecord(name, _literal(src, name)))
    if min_quality > 0:
        return filter_by_quality(records, min_quality)
    return records, []


def load_sequence(source: str, label: str = "sequence") -> SequenceRecord:
    """A single sequence (reference / query): first record of a file, or a literal."""
    if os.path.isfile(source):
        return read_records(source)[0]
    return SequenceRecord(label, _literal(source, label))


def format_fasta(records, width: int = 60) -> str:
    """records: iterable of SequenceRecord or (name, sequence) pairs."""
    out = []
    for r in records:
        name, seq = (r.name, r.sequence) if isinstance(r, SequenceRecord) else r
        out.append(f">{name}")
        out.extend(seq[i:i + width] for i in range(0, len(seq), width))
    return "\n".join(out) + "\n"


def write_text(path: str, text: str) -> None:
    directory = os.path.dirname(path)
    if directory:
        os.makedirs(directory, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)
