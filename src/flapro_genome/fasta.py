from __future__ import annotations

import gzip
import re
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

from .utils import FlaProGenomeError


AA_ALPHABET = set("ABCDEFGHIKLMNPQRSTVWXYZJUO*-?")
NUCLEOTIDE_ALPHABET = set("ACGTUNRYSWKMBDHV-?")
GENBANK_SUFFIXES = (".gb", ".gbf", ".gbk", ".gbff", ".genbank")


@dataclass(frozen=True)
class FastaRecord:
    id: str
    sequence: str
    description: str = ""


def _open_text(path: Path):
    return gzip.open(path, "rt") if path.name.endswith(".gz") else path.open("r")


def read_fasta(path: str | Path) -> Iterator[FastaRecord]:
    fasta_path = Path(path)
    try:
        handle = _open_text(fasta_path)
    except OSError as exc:
        raise FlaProGenomeError(f"Cannot open FASTA file {fasta_path}: {exc}") from exc
    with handle:
        header: str | None = None
        sequence: list[str] = []
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield _make_record(header, sequence, fasta_path)
                header = line[1:].strip()
                sequence = []
            elif header is None:
                raise FlaProGenomeError(
                    f"Invalid FASTA {fasta_path}: sequence before first header at line {line_number}"
                )
            else:
                sequence.append(re.sub(r"\s+", "", line).upper())
        if header is not None:
            yield _make_record(header, sequence, fasta_path)


def _make_record(header: str, sequence: list[str], path: Path) -> FastaRecord:
    if not header:
        raise FlaProGenomeError(f"Empty FASTA header in {path}")
    identifier = header.split()[0]
    seq = "".join(sequence)
    if not seq:
        raise FlaProGenomeError(f"Empty sequence for FASTA record {identifier!r} in {path}")
    return FastaRecord(identifier, seq, header)


def write_fasta(records: Iterable[FastaRecord], path: str | Path, width: int = 80) -> int:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with output.open("w") as handle:
        for record in records:
            handle.write(f">{record.description or record.id}\n")
            for start in range(0, len(record.sequence), width):
                handle.write(record.sequence[start : start + width] + "\n")
            count += 1
    return count


def validate_protein_sequence(sequence: str) -> bool:
    return bool(sequence) and set(sequence.upper()) <= AA_ALPHABET


def infer_sequence_file_format(path: str | Path) -> str:
    """Return ``fasta`` or ``genbank`` from a known suffix or file signature."""
    source = Path(path)
    name = source.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    if name.endswith(GENBANK_SUFFIXES):
        return "genbank"
    try:
        handle = _open_text(source)
    except OSError as exc:
        raise FlaProGenomeError(f"Cannot open sequence file {source}: {exc}") from exc
    with handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                return "fasta"
            if line.startswith("LOCUS"):
                return "genbank"
            break
    raise FlaProGenomeError(
        f"Could not identify {source} as FASTA or GenBank format"
    )


def infer_fasta_type(path: str | Path, max_residues: int = 100_000) -> str:
    observed: set[str] = set()
    n_residues = 0
    n_records = 0
    for record in read_fasta(path):
        n_records += 1
        observed.update(record.sequence.upper())
        n_residues += len(record.sequence)
        if n_residues >= max_residues:
            break
    if n_records == 0:
        raise FlaProGenomeError(f"No records found in FASTA file {path}")
    invalid = observed - AA_ALPHABET
    if invalid:
        raise FlaProGenomeError(
            f"FASTA {path} contains unsupported sequence characters: {''.join(sorted(invalid))}"
        )
    return "nucleotide" if observed <= NUCLEOTIDE_ALPHABET else "protein"
