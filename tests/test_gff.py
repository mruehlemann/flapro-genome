from __future__ import annotations

from pathlib import Path

import pytest

from flapro_genome.fasta import infer_sequence_file_format, read_fasta
from flapro_genome.gff import extract_gff
from flapro_genome.utils import FlaProGenomeError


def test_extract_gff_with_paired_fasta(tmp_path: Path) -> None:
    gff = tmp_path / "genome.gff3"
    fasta = tmp_path / "genome.fna"
    gff.write_text(
        "##gff-version 3\n"
        "contig_1\tProdigal\tCDS\t1\t12\t.\t+\t0\tID=cds1;protein_id=GFF_PROT_1\n"
    )
    fasta.write_text(">contig_1\nATGGCTGAATAA\n")

    result = extract_gff(gff, tmp_path / "proteins.faa", tmp_path / "contigs.fna", fasta)

    assert result.n_cds_features == 1
    assert result.n_translations_derived == 1
    assert [(record.id, record.sequence) for record in result.proteins] == [
        ("GFF_PROT_1", "MAE")
    ]


def test_extract_gff_with_embedded_fasta_and_reverse_strand(tmp_path: Path) -> None:
    gff = tmp_path / "genome.gff"
    gff.write_text(
        "##gff-version 3\n"
        "contig_1\ttool\tCDS\t1\t12\t.\t-\t0\tID=minus_cds;locus_tag=MINUS_1\n"
        "##FASTA\n"
        ">contig_1\n"
        "TTATTCAGCCAT\n"
    )

    result = extract_gff(gff, tmp_path / "proteins.faa", tmp_path / "contigs.fna")

    assert infer_sequence_file_format(gff) == "gff"
    assert result.proteins[0].id == "MINUS_1"
    assert result.proteins[0].sequence == "MAE"


def test_extract_gff_uses_translation_attribute_without_fasta(tmp_path: Path) -> None:
    gff = tmp_path / "translated.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "contig_1\ttool\tCDS\t1\t12\t.\t+\t0\tID=cds1;translation=MAE\n"
    )
    protein_path = tmp_path / "proteins.faa"
    result = extract_gff(gff, protein_path, tmp_path / "contigs.fna")
    assert result.n_translations_derived == 0
    assert result.contigs == []
    assert [(record.id, record.sequence) for record in read_fasta(protein_path)] == [
        ("cds1", "MAE")
    ]


def test_coordinate_only_gff_requires_sequence(tmp_path: Path) -> None:
    gff = tmp_path / "genome.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "contig_1\ttool\tCDS\t1\t12\t.\t+\t0\tID=cds1\n"
    )
    with pytest.raises(FlaProGenomeError, match="--sequence genome.fna"):
        extract_gff(gff, tmp_path / "proteins.faa", tmp_path / "contigs.fna")


def test_gff_phase_is_applied_before_translation(tmp_path: Path) -> None:
    gff = tmp_path / "phased.gff3"
    fasta = tmp_path / "phased.fna"
    gff.write_text(
        "##gff-version 3\n"
        "contig_1\ttool\tCDS\t1\t13\t.\t+\t1\tID=phased_cds\n"
    )
    fasta.write_text(">contig_1\nAATGGCTGAATAA\n")
    result = extract_gff(gff, tmp_path / "proteins.faa", tmp_path / "contigs.fna", fasta)
    assert result.proteins[0].sequence == "MAE"
