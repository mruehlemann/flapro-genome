from __future__ import annotations

import gzip
from pathlib import Path

from flapro_genome.fasta import infer_fasta_type, read_fasta


def test_infer_protein_and_preserve_id(tmp_path: Path) -> None:
    path = tmp_path / "proteins.faa"
    path.write_text(">protein_1 description here\nMPEPTIDEK\n")
    assert infer_fasta_type(path) == "protein"
    record = next(read_fasta(path))
    assert record.id == "protein_1"
    assert record.description == "protein_1 description here"


def test_infer_gzipped_nucleotide(tmp_path: Path) -> None:
    path = tmp_path / "genome.fna.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(">contig_1\nACGTNNRYACGT\n")
    assert infer_fasta_type(path) == "nucleotide"

