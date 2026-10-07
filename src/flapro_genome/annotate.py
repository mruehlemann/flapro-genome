from __future__ import annotations

import logging
import tempfile
from collections import Counter
from pathlib import Path

import pandas as pd

from .build_db import validate_database
from .classify import Thresholds, classify_hits
from .fasta import (
    FastaRecord,
    infer_fasta_type,
    infer_sequence_file_format,
    read_fasta,
    write_fasta,
)
from .genbank import extract_genbank
from .utils import FlaProGenomeError, require_executable, run_command, safe_stem


LOG = logging.getLogger("flapro_genome")
MMSEQS_COLUMNS = [
    "query", "target", "pident", "alnlen", "mismatch", "gapopen", "qstart", "qend",
    "qlen", "tstart", "tend", "tlen", "evalue", "bitscore",
]
MMSEQS_FORMAT = [*MMSEQS_COLUMNS[:-1], "bits"]
NUMERIC_COLUMNS = [
    "pident", "alnlen", "mismatch", "gapopen", "qstart", "qend", "qlen", "tstart",
    "tend", "tlen", "evalue", "bitscore",
]
ANNOTATION_COLUMNS = [
    "genome", "query_protein", "query_length", "classification",
    "nearest_flapro_protein", "pident", "alignment_length", "qcov", "tcov", "evalue", "bitscore",
    "assigned_cluster", "nearest_cluster", "cluster_consistent", "n_near_best_hits",
    "n_near_best_clusters", "near_best_clusters", "Cluster_CDHit",
    "reference_phylum", "reference_class", "reference_order", "reference_family",
    "reference_genus", "reference_species", "tlr5_class", "tlr5_source",
    "nearest_reference_tlr5_class", "cluster_phenotype_consistent", "tlr5_reference_protein",
    "tlr5_reference_identity", "reference_experimental_phenotype", "reference_marker_type",
]
RAW_HIT_COLUMNS = [
    "query", "target", "pident", "alnlen", "qlen", "tlen", "qcov", "tcov",
    "evalue", "bitscore", "target_cluster", "target_species", "target_predicted_v3",
]
SUMMARY_COLUMNS = [
    "genome", "n_proteins", "n_flagellin_candidates", "n_exact_reference",
    "n_known_flapro_type", "n_related_flagellin", "n_ambiguous", "n_weak_candidate",
    "n_unique_assigned_clusters", "n_tlr5_active", "n_tlr5_silent", "n_tlr5_unresolved",
]


def load_metadata(database_dir: str | Path) -> pd.DataFrame:
    path = Path(database_dir) / "flapro_reference_metadata.tsv"
    if not path.is_file():
        raise FlaProGenomeError(f"Reference metadata not found: {path}")
    metadata = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if metadata["Protein_ID"].duplicated().any():
        raise FlaProGenomeError(f"Reference metadata contains duplicate Protein_ID values: {path}")
    return metadata


def _normalize_fasta(input_path: Path, output_path: Path) -> list[FastaRecord]:
    records = list(read_fasta(input_path))
    ids = [record.id for record in records]
    if len(ids) != len(set(ids)):
        duplicates = sorted(identifier for identifier, count in Counter(ids).items() if count > 1)
        raise FlaProGenomeError("Input FASTA contains duplicate IDs: " + ", ".join(duplicates[:10]))
    write_fasta(records, output_path)
    return records


def predict_proteins(
    genome_fasta: Path,
    protein_output: Path,
    genes_output: Path,
    gff_output: Path,
    mode: str,
) -> None:
    require_executable("prodigal")
    run_command(
        [
            "prodigal", "-i", genome_fasta, "-a", protein_output, "-d", genes_output,
            "-f", "gff", "-o", gff_output, "-p", mode,
        ]
    )


def run_mmseqs_search(
    query_fasta: Path,
    database_dir: Path,
    work_dir: Path,
    threads: int,
    sensitivity: float = 7.5,
    max_seqs: int = 20,
    search_evalue: float = 1e-5,
) -> pd.DataFrame:
    require_executable("mmseqs")
    work_dir.mkdir(parents=True, exist_ok=True)
    ref_db = database_dir / "flapro_mmseqs" / "refDB"
    if not ref_db.with_suffix(".dbtype").is_file():
        raise FlaProGenomeError(f"MMseqs2 reference database not found: {ref_db}")
    query_db = work_dir / "queryDB"
    result_db = work_dir / "resultDB"
    mmseqs_tmp = work_dir / "tmp"
    hits_path = work_dir / "hits.tsv"
    run_command(["mmseqs", "createdb", query_fasta, query_db])
    run_command(
        [
            "mmseqs", "search", query_db, ref_db, result_db, mmseqs_tmp,
            "--threads", str(threads), "-s", str(sensitivity), "-e", str(search_evalue),
            "--max-seqs", str(max_seqs), "-a",
        ]
    )
    run_command(
        [
            "mmseqs", "convertalis", query_db, ref_db, result_db, hits_path,
            "--format-output", ",".join(MMSEQS_FORMAT),
        ]
    )
    if not hits_path.exists() or hits_path.stat().st_size == 0:
        return pd.DataFrame(columns=[*MMSEQS_COLUMNS, "qcov", "tcov"])
    hits = pd.read_csv(hits_path, sep="\t", names=MMSEQS_COLUMNS)
    for column in NUMERIC_COLUMNS:
        hits[column] = pd.to_numeric(hits[column], errors="raise")
    hits["qcov"] = hits["alnlen"] / hits["qlen"]
    hits["tcov"] = hits["alnlen"] / hits["tlen"]
    return hits


def enrich_raw_hits(hits: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    if hits.empty:
        return pd.DataFrame(columns=RAW_HIT_COLUMNS)
    lookup = metadata.set_index("Protein_ID")
    result = hits.copy()
    result["target_cluster"] = result["target"].map(lookup["Cluster_c4_representative"])
    result["target_species"] = result["target"].map(lookup["Species"])
    result["target_predicted_v3"] = result["target"].map(lookup["Predicted_v3"])
    return result[RAW_HIT_COLUMNS]


def make_summary(genome: str, n_proteins: int, annotations: pd.DataFrame) -> pd.DataFrame:
    if annotations.empty:
        counts: dict[str, int] = {}
        assigned_clusters = 0
        active = silent = unresolved = 0
    else:
        counts = annotations["classification"].value_counts().to_dict()
        assigned_clusters = int(annotations["assigned_cluster"].dropna().nunique())
        active = int((annotations["tlr5_class"] == "active").sum())
        silent = int((annotations["tlr5_class"] == "silent").sum())
        unresolved = int(len(annotations) - active - silent)
    row = {
        "genome": genome,
        "n_proteins": n_proteins,
        "n_flagellin_candidates": int(len(annotations)),
        "n_exact_reference": int(counts.get("exact_reference", 0)),
        "n_known_flapro_type": int(counts.get("known_flapro_type", 0)),
        "n_related_flagellin": int(counts.get("related_flagellin", 0)),
        "n_ambiguous": int(counts.get("ambiguous", 0)),
        "n_weak_candidate": int(counts.get("weak_candidate", 0)),
        "n_unique_assigned_clusters": assigned_clusters,
        "n_tlr5_active": active,
        "n_tlr5_silent": silent,
        "n_tlr5_unresolved": unresolved,
    }
    return pd.DataFrame([row], columns=SUMMARY_COLUMNS)


def annotate_fasta(
    input_fasta: str | Path,
    database_dir: str | Path,
    output_prefix: str | Path,
    thresholds: Thresholds | None = None,
    threads: int = 1,
    prodigal_mode: str = "meta",
    include_weak: bool = False,
) -> dict[str, Path]:
    thresholds = thresholds or Thresholds()
    source = Path(input_fasta).resolve()
    database = Path(database_dir).resolve()
    prefix = Path(output_prefix).resolve()
    if not source.is_file():
        raise FlaProGenomeError(f"Input sequence file not found: {source}")
    if threads < 1:
        raise FlaProGenomeError("--threads must be at least 1")
    errors = validate_database(database)
    if errors:
        raise FlaProGenomeError("Database validation failed:\n- " + "\n- ".join(errors))
    prefix.parent.mkdir(parents=True, exist_ok=True)
    genome = safe_stem(source)
    input_format = infer_sequence_file_format(source)

    with tempfile.TemporaryDirectory(prefix=f".{prefix.name}.work-", dir=prefix.parent) as temp_name:
        work = Path(temp_name)
        normalized_nucleotide = work / "input.fna"
        normalized_protein = work / "input.faa"
        use_prodigal = False
        if input_format == "genbank":
            extraction = extract_genbank(source, normalized_protein, normalized_nucleotide)
            if extraction.proteins:
                proteins = extraction.proteins
                protein_fasta = normalized_protein
                LOG.info("Using annotated CDS proteins from GenBank input %s", source)
            else:
                use_prodigal = True
                LOG.warning(
                    "No usable annotated CDS proteins found in %s; falling back to Prodigal",
                    source,
                )
        else:
            input_type = infer_fasta_type(source)
            LOG.info("Input %s inferred as %s FASTA", source, input_type)
            normalized_input = normalized_nucleotide if input_type == "nucleotide" else normalized_protein
            source_records = _normalize_fasta(source, normalized_input)
            if input_type == "nucleotide":
                use_prodigal = True
            else:
                proteins = source_records
                protein_fasta = normalized_protein

        if use_prodigal:
            protein_fasta = Path(str(prefix) + ".prodigal.faa")
            genes_fasta = Path(str(prefix) + ".prodigal.genes.fna")
            genes_gff = Path(str(prefix) + ".prodigal.gff")
            predict_proteins(
                normalized_nucleotide, protein_fasta, genes_fasta, genes_gff, prodigal_mode
            )
            proteins = list(read_fasta(protein_fasta))

        if proteins:
            hits = run_mmseqs_search(
                protein_fasta, database, work, threads, search_evalue=thresholds.candidate_evalue
            )
        else:
            LOG.warning("Prodigal predicted no proteins in %s", source)
            hits = pd.DataFrame(columns=[*MMSEQS_COLUMNS, "qcov", "tcov"])
        metadata = load_metadata(database)
        annotations = classify_hits(hits, metadata, thresholds, genome)
        if annotations.empty:
            annotations = pd.DataFrame(columns=ANNOTATION_COLUMNS)
        else:
            annotations = annotations.reindex(columns=ANNOTATION_COLUMNS)

        annotation_path = Path(str(prefix) + ".flagellins.tsv")
        raw_hits_path = Path(str(prefix) + ".all_hits.tsv")
        candidate_fasta_path = Path(str(prefix) + ".flagellins.faa")
        summary_path = Path(str(prefix) + ".summary.tsv")
        annotations.to_csv(annotation_path, sep="\t", index=False, na_rep="NA")
        enrich_raw_hits(hits, metadata).to_csv(raw_hits_path, sep="\t", index=False, na_rep="NA")

        selected_classes = {"exact_reference", "known_flapro_type", "related_flagellin"}
        if include_weak:
            selected_classes.add("weak_candidate")
        selected_ids = set(
            annotations.loc[annotations["classification"].isin(selected_classes), "query_protein"]
        )
        write_fasta((record for record in proteins if record.id in selected_ids), candidate_fasta_path)
        summary = make_summary(genome, len(proteins), annotations)
        summary.to_csv(summary_path, sep="\t", index=False, na_rep="NA")

    counts = summary.iloc[0]
    LOG.info("Detected %d candidate flagellins in %s", counts["n_flagellin_candidates"], genome)
    LOG.info(
        "%d exact references, %d known FlaPro types, %d related, %d ambiguous, %d weak",
        counts["n_exact_reference"], counts["n_known_flapro_type"],
        counts["n_related_flagellin"], counts["n_ambiguous"], counts["n_weak_candidate"],
    )
    return {
        "annotations": annotation_path,
        "raw_hits": raw_hits_path,
        "candidate_fasta": candidate_fasta_path,
        "summary": summary_path,
    }


def batch_annotate(
    input_dir: str | Path,
    database_dir: str | Path,
    output_dir: str | Path,
    thresholds: Thresholds | None = None,
    threads: int = 1,
    prodigal_mode: str = "meta",
    include_weak: bool = False,
) -> tuple[Path, Path]:
    source_dir = Path(input_dir).resolve()
    output = Path(output_dir).resolve()
    if not source_dir.is_dir():
        raise FlaProGenomeError(f"Batch input directory not found: {source_dir}")
    allowed = (
        ".faa", ".fasta", ".fa", ".fna", ".faa.gz", ".fasta.gz", ".fa.gz", ".fna.gz",
        ".gb", ".gbf", ".gbk", ".gbff", ".genbank", ".gb.gz", ".gbf.gz", ".gbk.gz",
        ".gbff.gz", ".genbank.gz",
    )
    inputs = sorted(path for path in source_dir.iterdir() if path.is_file() and path.name.lower().endswith(allowed))
    if not inputs:
        raise FlaProGenomeError(f"No supported FASTA or GenBank files found in {source_dir}")
    names = [safe_stem(path) for path in inputs]
    duplicate_names = sorted({name for name in names if names.count(name) > 1})
    if duplicate_names:
        raise FlaProGenomeError(
            "Batch inputs resolve to duplicate genome names: " + ", ".join(duplicate_names)
        )
    output.mkdir(parents=True, exist_ok=True)
    annotation_frames: list[pd.DataFrame] = []
    summary_frames: list[pd.DataFrame] = []
    for path in inputs:
        name = safe_stem(path)
        sample_dir = output / name
        sample_dir.mkdir(parents=True, exist_ok=True)
        paths = annotate_fasta(
            path, database_dir, sample_dir / name, thresholds, threads, prodigal_mode, include_weak
        )
        annotation_frames.append(pd.read_csv(paths["annotations"], sep="\t"))
        summary_frames.append(pd.read_csv(paths["summary"], sep="\t"))
    combined_annotations = output / "all_flagellins.tsv"
    combined_summaries = output / "all_genome_summaries.tsv"
    pd.concat(annotation_frames, ignore_index=True).to_csv(
        combined_annotations, sep="\t", index=False, na_rep="NA"
    )
    pd.concat(summary_frames, ignore_index=True).to_csv(
        combined_summaries, sep="\t", index=False, na_rep="NA"
    )
    return combined_annotations, combined_summaries
