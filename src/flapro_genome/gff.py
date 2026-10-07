from __future__ import annotations

import gzip
import logging
import re
from collections import Counter, OrderedDict
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from urllib.parse import unquote

from Bio import SeqIO
from Bio.Seq import Seq

from .fasta import (
    FastaRecord,
    infer_fasta_type,
    read_fasta,
    validate_protein_sequence,
    write_fasta,
)
from .utils import FlaProGenomeError


LOG = logging.getLogger("flapro_genome")


@dataclass(frozen=True)
class GFFPart:
    seqid: str
    start: int
    end: int
    strand: str
    phase: int
    attributes: dict[str, str]
    line_number: int


@dataclass(frozen=True)
class GFFExtraction:
    proteins: list[FastaRecord]
    contigs: list[FastaRecord]
    n_cds_features: int
    n_translations_derived: int


def _open_gff(path: Path):
    return gzip.open(path, "rt") if path.name.lower().endswith(".gz") else path.open("r")


def _parse_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for field in text.strip().strip(";").split(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field:
            key, value = field.split("=", 1)
        elif " " in field:  # tolerate common GTF-style key "value" fields
            key, value = field.split(None, 1)
        else:
            key, value = field, "true"
        attributes[unquote(key.strip())] = unquote(value.strip().strip('"'))
    return attributes


def _attribute(attributes: dict[str, str], *names: str) -> str | None:
    lower = {key.lower(): value for key, value in attributes.items()}
    for name in names:
        value = lower.get(name.lower())
        if value:
            return value.split(",", 1)[0].strip()
    return None


def _is_pseudogene(attributes: dict[str, str]) -> bool:
    lower = {key.lower(): value.lower() for key, value in attributes.items()}
    false_values = {"0", "false", "no"}
    pseudo = "pseudo" in lower and lower["pseudo"] not in false_values
    pseudogene = "pseudogene" in lower and lower["pseudogene"] not in false_values
    return pseudo or pseudogene or lower.get("gbkey") == "pseudogene"


def _group_key(part: GFFPart, ordinal: int) -> str:
    return (
        _attribute(part.attributes, "protein_id")
        or _attribute(part.attributes, "Parent")
        or _attribute(part.attributes, "ID")
        or f"{part.seqid}:{part.start}-{part.end}:{part.strand}:{ordinal}"
    )


def _protein_identifier(parts: list[GFFPart], ordinal: int) -> str:
    first = parts[0]
    identifier = (
        _attribute(first.attributes, "protein_id")
        or _attribute(first.attributes, "locus_tag")
        or _attribute(first.attributes, "Name")
        or _attribute(first.attributes, "gene")
        or (_attribute(first.attributes, "Parent") if len(parts) > 1 else None)
        or _attribute(first.attributes, "ID")
        or _attribute(first.attributes, "Parent")
        or f"{first.seqid}_CDS_{ordinal}"
    )
    sanitized = re.sub(r"\s+", "_", identifier)
    if sanitized != identifier:
        LOG.warning("Replaced whitespace in GFF protein identifier %r", identifier)
    return sanitized


def _read_gff(path: Path) -> tuple[list[GFFPart], list[FastaRecord]]:
    parts: list[GFFPart] = []
    fasta_lines: list[str] = []
    in_fasta = False
    try:
        with _open_gff(path) as handle:
            for line_number, raw in enumerate(handle, start=1):
                if in_fasta:
                    fasta_lines.append(raw)
                    continue
                line = raw.rstrip("\n\r")
                if line == "##FASTA":
                    in_fasta = True
                    continue
                if not line or line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) != 9:
                    raise FlaProGenomeError(
                        f"Invalid GFF row at {path}:{line_number}: expected 9 tab-separated fields"
                    )
                seqid, _source, feature_type, start, end, _score, strand, phase, attribute_text = fields
                if feature_type.lower() != "cds":
                    continue
                try:
                    start_value = int(start)
                    end_value = int(end)
                    phase_value = 0 if phase == "." else int(phase)
                except ValueError as exc:
                    raise FlaProGenomeError(
                        f"Invalid CDS coordinates or phase at {path}:{line_number}"
                    ) from exc
                if start_value < 1 or end_value < start_value:
                    raise FlaProGenomeError(f"Invalid CDS interval at {path}:{line_number}")
                if strand not in {"+", "-"}:
                    raise FlaProGenomeError(
                        f"CDS strand must be '+' or '-' at {path}:{line_number}"
                    )
                if phase_value not in {0, 1, 2}:
                    raise FlaProGenomeError(f"CDS phase must be 0, 1, or 2 at {path}:{line_number}")
                attributes = _parse_attributes(attribute_text)
                if _is_pseudogene(attributes):
                    continue
                parts.append(GFFPart(
                    unquote(seqid), start_value, end_value, strand, phase_value,
                    attributes, line_number,
                ))
    except OSError as exc:
        raise FlaProGenomeError(f"Cannot open GFF file {path}: {exc}") from exc

    embedded: list[FastaRecord] = []
    if fasta_lines:
        try:
            embedded = [
                FastaRecord(str(record.id), str(record.seq).upper(), str(record.description))
                for record in SeqIO.parse(StringIO("".join(fasta_lines)), "fasta")
            ]
        except Exception as exc:
            raise FlaProGenomeError(f"Cannot parse embedded FASTA section in {path}: {exc}") from exc
        if not embedded:
            raise FlaProGenomeError(f"GFF file {path} has an empty ##FASTA section")
    return parts, embedded


def _load_external_sequences(path: Path) -> list[FastaRecord]:
    if infer_fasta_type(path) != "nucleotide":
        raise FlaProGenomeError(f"Paired GFF sequence file must be nucleotide FASTA: {path}")
    return list(read_fasta(path))


def _derive_translation(parts: list[GFFPart], sequences: dict[str, str]) -> str:
    seqids = {part.seqid for part in parts}
    strands = {part.strand for part in parts}
    if len(seqids) != 1 or len(strands) != 1:
        raise FlaProGenomeError("Grouped GFF CDS rows span multiple sequences or strands")
    seqid = parts[0].seqid
    if seqid not in sequences:
        raise FlaProGenomeError(
            f"GFF sequence ID {seqid!r} is absent from the embedded or paired FASTA"
        )
    genome = sequences[seqid]
    ordered = sorted(parts, key=lambda part: part.start, reverse=parts[0].strand == "-")
    fragments: list[str] = []
    for part in ordered:
        if part.end > len(genome):
            raise FlaProGenomeError(
                f"GFF CDS at line {part.line_number} ends beyond sequence {seqid!r}"
            )
        fragment = Seq(genome[part.start - 1 : part.end])
        if part.strand == "-":
            fragment = fragment.reverse_complement()
        fragments.append(str(fragment))
    coding = "".join(fragments)[ordered[0].phase :]
    remainder = len(coding) % 3
    if remainder:
        coding = coding[:-remainder]
    if not coding:
        raise FlaProGenomeError(f"GFF CDS on {seqid!r} has no complete codon")
    table_text = _attribute(ordered[0].attributes, "transl_table", "translation_table") or "11"
    try:
        table = int(table_text)
        return str(Seq(coding).translate(table=table, to_stop=False)).rstrip("*")
    except Exception as exc:
        raise FlaProGenomeError(
            f"Could not translate GFF CDS on {seqid!r}: {exc}"
        ) from exc


def extract_gff(
    input_path: str | Path,
    protein_output: str | Path,
    nucleotide_output: str | Path,
    sequence_path: str | Path | None = None,
) -> GFFExtraction:
    """Extract CDS proteins from GFF3 plus embedded or paired nucleotide sequences."""
    source = Path(input_path)
    parts, embedded_contigs = _read_gff(source)
    if sequence_path is not None and embedded_contigs:
        LOG.info("Using --sequence FASTA instead of the embedded GFF FASTA section")
    contigs = (
        _load_external_sequences(Path(sequence_path))
        if sequence_path is not None
        else embedded_contigs
    )
    contig_ids = [record.id for record in contigs]
    if len(contig_ids) != len(set(contig_ids)):
        raise FlaProGenomeError("GFF nucleotide sequences contain duplicate IDs")
    sequences = {record.id: record.sequence for record in contigs}

    groups: OrderedDict[str, list[GFFPart]] = OrderedDict()
    for ordinal, part in enumerate(parts, start=1):
        groups.setdefault(_group_key(part, ordinal), []).append(part)

    proteins: list[FastaRecord] = []
    n_derived = 0
    missing_sequence_groups: list[str] = []
    for ordinal, group_parts in enumerate(groups.values(), start=1):
        identifier = _protein_identifier(group_parts, ordinal)
        translations = {
            re.sub(r"\s+", "", value).upper().rstrip("*")
            for part in group_parts
            if (value := _attribute(
                part.attributes, "translation", "protein_sequence", "aa_sequence"
            ))
        }
        if len(translations) > 1:
            raise FlaProGenomeError(f"Conflicting GFF translation attributes for CDS {identifier}")
        if translations:
            translation = translations.pop()
        elif sequences:
            translation = _derive_translation(group_parts, sequences)
            n_derived += 1
        else:
            missing_sequence_groups.append(identifier)
            continue
        if not translation:
            continue
        if not validate_protein_sequence(translation):
            raise FlaProGenomeError(f"Invalid amino-acid characters in GFF CDS {identifier}")
        proteins.append(FastaRecord(identifier, translation))

    if missing_sequence_groups:
        raise FlaProGenomeError(
            f"GFF file {source} has CDS coordinates without translations or nucleotide sequence. "
            "Provide the matching genome with --sequence genome.fna. "
            f"First affected CDS: {missing_sequence_groups[0]}"
        )
    protein_ids = [record.id for record in proteins]
    if len(protein_ids) != len(set(protein_ids)):
        duplicates = sorted(
            identifier for identifier, count in Counter(protein_ids).items() if count > 1
        )
        raise FlaProGenomeError(
            "GFF CDS features resolve to duplicate protein identifiers: " + ", ".join(duplicates[:10])
        )

    write_fasta(proteins, protein_output)
    write_fasta(contigs, nucleotide_output)
    LOG.info(
        "Read %d GFF CDS row(s), produced %d protein(s), and derived %d translation(s)",
        len(parts), len(proteins), n_derived,
    )
    return GFFExtraction(proteins, contigs, len(parts), n_derived)
