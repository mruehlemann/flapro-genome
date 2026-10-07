from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from flapro_genome.annotate import MMSEQS_COLUMNS, annotate_fasta
from flapro_genome.fasta import FastaRecord, infer_sequence_file_format, read_fasta, write_fasta
from flapro_genome.genbank import extract_genbank


def _write_genbank(path: Path, include_translation: bool = True) -> None:
    sequence = Seq("ATGGCTGAATAA")
    record = SeqRecord(sequence, id="contig_1", name="contig_1", description="test record")
    record.annotations["molecule_type"] = "DNA"
    qualifiers = {
        "locus_tag": ["TEST_0001"],
        "protein_id": ["WP_TEST_1"],
        "transl_table": ["11"],
    }
    if include_translation:
        qualifiers["translation"] = ["MAE"]
    record.features.append(
        SeqFeature(FeatureLocation(0, len(sequence), strand=1), type="CDS", qualifiers=qualifiers)
    )
    SeqIO.write([record], path, "genbank")


def test_extract_genbank_uses_annotated_translation_and_protein_id(tmp_path: Path) -> None:
    source = tmp_path / "genome.gbk"
    proteins = tmp_path / "proteins.faa"
    nucleotides = tmp_path / "contigs.fna"
    _write_genbank(source)

    result = extract_genbank(source, proteins, nucleotides)

    assert result.n_cds_features == 1
    assert result.n_translations_derived == 0
    assert [(record.id, record.sequence) for record in result.proteins] == [("WP_TEST_1", "MAE")]
    assert [(record.id, record.sequence) for record in read_fasta(nucleotides)] == [
        ("contig_1", "ATGGCTGAATAA")
    ]


def test_extract_genbank_derives_missing_translation(tmp_path: Path) -> None:
    source = tmp_path / "genome.gbff"
    _write_genbank(source, include_translation=False)
    result = extract_genbank(source, tmp_path / "proteins.faa", tmp_path / "contigs.fna")
    assert result.n_translations_derived == 1
    assert result.proteins[0].sequence == "MAE"


def test_detect_gzipped_genbank(tmp_path: Path) -> None:
    source = tmp_path / "genome.gbk"
    compressed = tmp_path / "genome.gbk.gz"
    _write_genbank(source)
    with source.open("rb") as input_handle, gzip.open(compressed, "wb") as output_handle:
        output_handle.write(input_handle.read())
    assert infer_sequence_file_format(compressed) == "genbank"
    result = extract_genbank(
        compressed, tmp_path / "compressed.faa", tmp_path / "compressed.fna"
    )
    assert result.proteins[0].id == "WP_TEST_1"


def test_annotate_genbank_uses_cds_proteins(tmp_path: Path, monkeypatch) -> None:
    source = tmp_path / "genome.gbk"
    _write_genbank(source)
    metadata = pd.DataFrame(
        [{
            "Protein_ID": "REF_A",
            "Cluster_c4_representative": "REF_A",
            "Cluster_CDHit": "1",
            "Phylum": "P",
            "Class": "C",
            "Order": "O",
            "Family": "F",
            "Genus": "G",
            "Species": "S",
            "Predicted_v3": "active",
            "Experimental": "not_checked",
            "MarkerType": "True Marker",
            "cluster_phenotype_consistent": "True",
        }]
    )

    monkeypatch.setattr("flapro_genome.annotate.validate_database", lambda _database: [])
    monkeypatch.setattr("flapro_genome.annotate.load_metadata", lambda _database: metadata)

    def fake_search(query_fasta, _database, _work, _threads, **_kwargs):
        records = list(read_fasta(query_fasta))
        assert [(record.id, record.sequence) for record in records] == [("WP_TEST_1", "MAE")]
        return pd.DataFrame(
            [{
                "query": "WP_TEST_1", "target": "REF_A", "pident": 100.0,
                "alnlen": 3, "mismatch": 0, "gapopen": 0, "qstart": 1,
                "qend": 3, "qlen": 3, "tstart": 1, "tend": 3, "tlen": 3,
                "evalue": 1e-20, "bitscore": 100.0, "qcov": 1.0, "tcov": 1.0,
            }]
        )

    monkeypatch.setattr("flapro_genome.annotate.run_mmseqs_search", fake_search)
    paths = annotate_fasta(source, tmp_path / "db", tmp_path / "result")
    annotations = pd.read_csv(paths["annotations"], sep="\t")
    assert annotations.loc[0, "query_protein"] == "WP_TEST_1"
    assert annotations.loc[0, "classification"] == "exact_reference"
    assert annotations.loc[0, "tlr5_class"] == "active"


def test_genbank_without_cds_falls_back_to_prodigal(tmp_path: Path, monkeypatch) -> None:
    source = tmp_path / "unannotated.gbk"
    record = SeqRecord(Seq("ATGGCTGAATAA"), id="contig_1", description="unannotated")
    record.annotations["molecule_type"] = "DNA"
    SeqIO.write([record], source, "genbank")

    monkeypatch.setattr("flapro_genome.annotate.validate_database", lambda _database: [])
    monkeypatch.setattr(
        "flapro_genome.annotate.load_metadata",
        lambda _database: pd.DataFrame(columns=["Protein_ID"]),
    )

    def fake_prodigal(_genome, protein_output, _genes_output, _gff_output, _mode):
        write_fasta([FastaRecord("predicted_1", "MAE")], protein_output)

    monkeypatch.setattr("flapro_genome.annotate.predict_proteins", fake_prodigal)
    monkeypatch.setattr(
        "flapro_genome.annotate.run_mmseqs_search",
        lambda *_args, **_kwargs: pd.DataFrame(columns=[*MMSEQS_COLUMNS, "qcov", "tcov"]),
    )
    paths = annotate_fasta(source, tmp_path / "db", tmp_path / "fallback")
    summary = pd.read_csv(paths["summary"], sep="\t")
    assert summary.loc[0, "n_proteins"] == 1
    assert summary.loc[0, "n_flagellin_candidates"] == 0
    assert (tmp_path / "fallback.prodigal.faa").is_file()
