from __future__ import annotations

import gzip
import logging
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature

from .fasta import AA_ALPHABET, FastaRecord, validate_protein_sequence, write_fasta
from .utils import FlaProGenomeError


LOG = logging.getLogger("flapro_genome")


@dataclass(frozen=True)
class GenBankExtraction:
    proteins: list[FastaRecord]
    contigs: list[FastaRecord]
    n_cds_features: int
    n_translations_derived: int


def _open_genbank(path: Path):
    return gzip.open(path, "rt") if path.name.lower().endswith(".gz") else path.open("r")


def _qualifier(feature: SeqFeature, name: str) -> str | None:
    values = feature.qualifiers.get(name, [])
    if not values:
        return None
    value = str(values[0]).strip()
    return value or None


def _protein_identifier(feature: SeqFeature, record_id: str, cds_number: int) -> str:
    identifier = (
        _qualifier(feature, "protein_id")
        or _qualifier(feature, "locus_tag")
        or _qualifier(feature, "gene")
        or f"{record_id}_CDS_{cds_number}"
    )
    sanitized = re.sub(r"\s+", "_", identifier)
    if sanitized != identifier:
        LOG.warning("Replaced whitespace in GenBank protein identifier %r", identifier)
    return sanitized


def _derive_translation(feature: SeqFeature, nucleotide_sequence) -> str | None:
    try:
        coding = feature.extract(nucleotide_sequence)
        codon_start = int(_qualifier(feature, "codon_start") or "1")
        if codon_start not in {1, 2, 3}:
            raise ValueError(f"invalid codon_start={codon_start}")
        coding = coding[codon_start - 1 :]
        remainder = len(coding) % 3
        if remainder:
            coding = coding[:-remainder]
        if not coding:
            return None
        table = int(_qualifier(feature, "transl_table") or "11")
        return str(coding.translate(table=table, to_stop=False)).rstrip("*")
    except Exception as exc:
        LOG.warning("Could not translate CDS feature at %s: %s", feature.location, exc)
        return None


def extract_genbank(
    input_path: str | Path,
    protein_output: str | Path,
    nucleotide_output: str | Path,
) -> GenBankExtraction:
    """Extract annotated CDS proteins and nucleotide records from GenBank."""
    source = Path(input_path)
    try:
        with _open_genbank(source) as handle:
            records = list(SeqIO.parse(handle, "genbank"))
    except Exception as exc:
        raise FlaProGenomeError(f"Cannot parse GenBank file {source}: {exc}") from exc
    if not records:
        raise FlaProGenomeError(f"No records found in GenBank file {source}")

    contigs: list[FastaRecord] = []
    proteins: list[FastaRecord] = []
    n_cds = 0
    n_derived = 0
    for record in records:
        contigs.append(FastaRecord(str(record.id), str(record.seq).upper()))
        for feature in record.features:
            if feature.type != "CDS":
                continue
            n_cds += 1
            if "pseudo" in feature.qualifiers or "pseudogene" in feature.qualifiers:
                continue
            identifier = _protein_identifier(feature, str(record.id), n_cds)
            translation = _qualifier(feature, "translation")
            if translation is not None:
                translation = re.sub(r"\s+", "", translation).upper().rstrip("*")
            else:
                translation = _derive_translation(feature, record.seq)
                if translation:
                    n_derived += 1
            if not translation:
                continue
            if not validate_protein_sequence(translation):
                invalid = "".join(sorted(set(translation) - AA_ALPHABET))
                raise FlaProGenomeError(
                    f"Invalid amino-acid characters ({invalid}) in GenBank CDS {identifier}"
                )
            proteins.append(FastaRecord(identifier, translation))

    contig_ids = [record.id for record in contigs]
    if len(contig_ids) != len(set(contig_ids)):
        raise FlaProGenomeError("GenBank file contains duplicate record IDs")
    protein_ids = [record.id for record in proteins]
    if len(protein_ids) != len(set(protein_ids)):
        duplicates = sorted(
            identifier for identifier, count in Counter(protein_ids).items() if count > 1
        )
        raise FlaProGenomeError(
            "GenBank CDS features contain duplicate protein identifiers: " + ", ".join(duplicates[:10])
        )

    write_fasta(contigs, nucleotide_output)
    write_fasta(proteins, protein_output)
    LOG.info(
        "Read %d GenBank record(s), %d CDS feature(s), and %d usable protein(s) (%d translations derived)",
        len(records), n_cds, len(proteins), n_derived,
    )
    return GenBankExtraction(proteins, contigs, n_cds, n_derived)
