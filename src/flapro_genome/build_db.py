from __future__ import annotations

import csv
import json
import logging
import statistics
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from . import __version__
from .fasta import FastaRecord, read_fasta, validate_protein_sequence, write_fasta
from .utils import FlaProGenomeError, command_version, require_executable, run_command


LOG = logging.getLogger("flapro_genome")
PARAMETERS_RELATIVE = Path("data/model_in/parameters_for_model_curated_2569_no_na_31-01-25.tsv")
TAXONOMY_RELATIVE = Path("data/taxonomy_cluster_repr_c4-pred3.tsv")
ID_EQUIVALENTS = {"proteinid", "proteinidentifier", "flagellinid", "id", "accession"}
METADATA_COLUMNS = [
    "Protein_ID",
    "length",
    "Cluster_CDHit",
    "Cluster_c4_representative",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species",
    "Predicted_v1",
    "Experimental",
    "Predicted_v3",
    "MarkerType",
    "Content_of_the_cluster",
    "cluster_n_active",
    "cluster_n_silent",
    "cluster_n_total",
    "cluster_phenotype_consistent",
    "cluster_representative_in_reference",
]


def _normalized_column(name: str) -> str:
    return "".join(character for character in name.lower() if character.isalnum())


def parse_parameter_tsv(path: str | Path) -> pd.DataFrame:
    """Parse the known FlaPro off-by-one parameter table safely."""
    source = Path(path)
    with source.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = [field.strip() for field in next(reader)]
        except StopIteration as exc:
            raise FlaProGenomeError(f"Parameter TSV is empty: {source}") from exc
        rows = [[field.strip() for field in row] for row in reader if any(field.strip() for field in row)]

    if not rows:
        raise FlaProGenomeError(f"Parameter TSV contains a header but no data: {source}")
    first_count = len(rows[0])
    header_count = len(header)
    if first_count == header_count + 1:
        header = ["Protein_ID", *header]
    elif first_count == header_count:
        candidates = [name for name in header if _normalized_column(name) in ID_EQUIVALENTS]
        if len(candidates) != 1:
            raise FlaProGenomeError(
                f"Parameter TSV has {header_count} header and data fields, but no unique "
                "Protein_ID-equivalent column. Columns: " + ", ".join(header)
            )
        header[header.index(candidates[0])] = "Protein_ID"
    else:
        raise FlaProGenomeError(
            f"Unexpected parameter TSV layout: header has {header_count} fields but first "
            f"data row has {first_count}; expected equal counts or exactly one extra data field"
        )
    if len(header) != len(set(header)):
        duplicates = sorted({name for name in header if header.count(name) > 1})
        raise FlaProGenomeError(
            "Parameter TSV resolves to duplicate column names: " + ", ".join(duplicates)
        )

    bad_rows = [index + 2 for index, row in enumerate(rows) if len(row) != len(header)]
    if bad_rows:
        preview = ", ".join(map(str, bad_rows[:10]))
        raise FlaProGenomeError(
            f"Parameter TSV has rows with unexpected field counts at line(s): {preview}"
        )

    frame = pd.DataFrame(rows, columns=header)
    required = {"Protein_ID", "Cluster_CDHit", "sequence"}
    missing = required - set(frame.columns)
    if missing:
        raise FlaProGenomeError(f"Parameter TSV is missing required columns: {', '.join(sorted(missing))}")
    if frame["Protein_ID"].eq("").any():
        raise FlaProGenomeError("Parameter TSV contains an empty Protein_ID")
    duplicates = frame.loc[frame["Protein_ID"].duplicated(keep=False), "Protein_ID"].unique()
    if len(duplicates):
        raise FlaProGenomeError(f"Duplicate Protein_ID values in parameter TSV: {', '.join(duplicates[:10])}")
    if frame["sequence"].eq("").any():
        empty_ids = frame.loc[frame["sequence"].eq(""), "Protein_ID"].head(10).tolist()
        raise FlaProGenomeError(f"Empty protein sequence(s) for: {', '.join(empty_ids)}")
    invalid_ids = frame.loc[~frame["sequence"].map(validate_protein_sequence), "Protein_ID"].head(10).tolist()
    if invalid_ids:
        raise FlaProGenomeError(f"Invalid amino-acid sequence(s) for: {', '.join(invalid_ids)}")
    frame["length"] = frame["sequence"].str.len()
    lengths = frame["length"].tolist()
    LOG.info(
        "Loaded %d full-length model sequences (length min=%d, median=%.1f, max=%d aa)",
        len(frame), min(lengths), statistics.median(lengths), max(lengths),
    )
    if not (1_500 <= len(frame) <= 4_000):
        LOG.warning("Parameter sequence count (%d) differs substantially from the expected ~2462", len(frame))
    if not (150 <= min(lengths) <= 300):
        LOG.warning("Minimum sequence length (%d aa) differs from the expected ~201 aa", min(lengths))
    if not (350 <= statistics.median(lengths) <= 650):
        LOG.warning(
            "Median sequence length (%.1f aa) differs from the expected ~469 aa",
            statistics.median(lengths),
        )
    if not (1_000 <= max(lengths) <= 2_500):
        LOG.warning("Maximum sequence length (%d aa) differs from the expected ~1623 aa", max(lengths))
    return frame


def load_taxonomy_tsv(path: str | Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {
        "Flagellin_ID", "Phylum", "Class", "Order", "Family", "Genus", "Species",
        "Predicted_v1", "Experimental", "Predicted_v3", "Cluster_c4_representative",
        "MarkerType", "Content_of_the_cluster",
    }
    missing = required - set(frame.columns)
    if missing:
        raise FlaProGenomeError(f"Taxonomy TSV is missing required columns: {', '.join(sorted(missing))}")
    duplicates = frame.loc[frame["Flagellin_ID"].duplicated(keep=False), "Flagellin_ID"].unique()
    if len(duplicates):
        raise FlaProGenomeError(f"Duplicate Flagellin_ID values in taxonomy TSV: {', '.join(duplicates[:10])}")
    return frame


def join_references(parameters: pd.DataFrame, taxonomy: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    parameter_ids = set(parameters["Protein_ID"])
    taxonomy_ids = set(taxonomy["Flagellin_ID"])
    shared = parameter_ids & taxonomy_ids
    id_stats = {
        "n_ids_only_parameter": len(parameter_ids - taxonomy_ids),
        "n_ids_only_taxonomy": len(taxonomy_ids - parameter_ids),
        "n_ids_both": len(shared),
    }
    LOG.info(
        "ID overlap: %d in both, %d only in parameter TSV, %d only in taxonomy TSV",
        id_stats["n_ids_both"], id_stats["n_ids_only_parameter"], id_stats["n_ids_only_taxonomy"],
    )
    joined = parameters.merge(taxonomy, left_on="Protein_ID", right_on="Flagellin_ID", how="inner", validate="one_to_one")
    joined["Predicted_v3"] = joined["Predicted_v3"].str.lower()
    excluded = joined.loc[~joined["Predicted_v3"].isin(["active", "silent"])]
    if not excluded.empty:
        LOG.warning(
            "Excluded %d joined proteins without an active/silent Predicted_v3 phenotype",
            len(excluded),
        )
    joined = joined.loc[joined["Predicted_v3"].isin(["active", "silent"])].copy()
    if len(parameters) >= 1_000 and not (1_500 <= len(joined) <= 4_000):
        LOG.warning(
            "Classified reference count (%d) differs substantially from the expected ~2439",
            len(joined),
        )

    group = joined.groupby("Cluster_c4_representative", dropna=False)["Predicted_v3"]
    cluster_stats = group.agg(
        cluster_n_active=lambda values: int((values == "active").sum()),
        cluster_n_silent=lambda values: int((values == "silent").sum()),
        cluster_n_total="size",
    ).reset_index()
    cluster_stats["cluster_phenotype_consistent"] = (
        (cluster_stats["cluster_n_active"] == 0) | (cluster_stats["cluster_n_silent"] == 0)
    )
    joined = joined.merge(cluster_stats, on="Cluster_c4_representative", how="left", validate="many_to_one")
    reference_ids = set(joined["Protein_ID"])
    joined["cluster_representative_in_reference"] = joined["Cluster_c4_representative"].isin(reference_ids)
    return joined, id_stats


def build_database(flapro_repo: str | Path, output_dir: str | Path) -> dict[str, object]:
    repo = Path(flapro_repo).resolve()
    output = Path(output_dir).resolve()
    parameter_path = repo / PARAMETERS_RELATIVE
    taxonomy_path = repo / TAXONOMY_RELATIVE
    missing = [path for path in (parameter_path, taxonomy_path) if not path.is_file()]
    if missing:
        raise FlaProGenomeError("Required FlaPro source file(s) missing: " + ", ".join(map(str, missing)))
    require_executable("mmseqs")

    parameters = parse_parameter_tsv(parameter_path)
    taxonomy = load_taxonomy_tsv(taxonomy_path)
    references, id_stats = join_references(parameters, taxonomy)
    output.mkdir(parents=True, exist_ok=True)
    fasta_path = output / "flapro_reference.faa"
    metadata_path = output / "flapro_reference_metadata.tsv"
    mmseqs_prefix = output / "flapro_mmseqs" / "refDB"
    mmseqs_prefix.parent.mkdir(parents=True, exist_ok=True)

    write_fasta(
        (FastaRecord(row.Protein_ID, row.sequence) for row in references.itertuples(index=False)),
        fasta_path,
    )
    references[METADATA_COLUMNS].to_csv(metadata_path, sep="\t", index=False, na_rep="NA")
    run_command(["mmseqs", "createdb", fasta_path, mmseqs_prefix])

    lengths = parameters["length"].astype(int).tolist()
    stats: dict[str, object] = {
        "n_parameter_sequences": int(len(parameters)),
        "n_taxonomy_entries": int(len(taxonomy)),
        "n_joined_references": int(len(references)),
        "n_active": int((references["Predicted_v3"] == "active").sum()),
        "n_silent": int((references["Predicted_v3"] == "silent").sum()),
        "n_clusters": int(references["Cluster_c4_representative"].nunique()),
        "n_mixed_phenotype_clusters": int(
            (~references.drop_duplicates("Cluster_c4_representative")["cluster_phenotype_consistent"]).sum()
        ),
        "sequence_length_min": int(min(lengths)),
        "sequence_length_median": float(statistics.median(lengths)),
        "sequence_length_max": int(max(lengths)),
        **id_stats,
    }
    (output / "build_stats.json").write_text(json.dumps(stats, indent=2) + "\n")

    commit = "unknown"
    try:
        commit = run_command(["git", "-C", repo, "rev-parse", "HEAD"]).stdout.strip() or "unknown"
    except FlaProGenomeError:
        pass
    version_text = "\n".join(
        [
            f"FlaPro_commit={commit}",
            f"database_build_date={datetime.now(timezone.utc).isoformat()}",
            f"flapro_genome_version={__version__}",
            f"MMseqs2_version={command_version(['mmseqs', 'version'])}",
        ]
    )
    (output / "VERSION").write_text(version_text + "\n")
    LOG.info("Matched %d sequences to classified FlaPro references", len(references))
    LOG.info("Built database containing %d proteins at %s", len(references), output)
    return stats


def validate_database(database_dir: str | Path) -> list[str]:
    database = Path(database_dir)
    errors: list[str] = []
    fasta_path = database / "flapro_reference.faa"
    metadata_path = database / "flapro_reference_metadata.tsv"
    for required in (fasta_path, metadata_path, database / "build_stats.json", database / "VERSION"):
        if not required.is_file():
            errors.append(f"Missing required database file: {required}")
    if errors:
        return errors
    try:
        records = list(read_fasta(fasta_path))
    except FlaProGenomeError as exc:
        return [str(exc)]
    fasta_ids = [record.id for record in records]
    if not fasta_ids:
        errors.append("Reference FASTA contains no sequences")
    if len(fasta_ids) != len(set(fasta_ids)):
        errors.append("Reference FASTA contains duplicate IDs")
    bad = [record.id for record in records if not validate_protein_sequence(record.sequence)]
    if bad:
        errors.append("Reference FASTA contains invalid/empty sequences: " + ", ".join(bad[:10]))
    try:
        metadata = pd.read_csv(metadata_path, sep="\t", dtype=str, keep_default_na=False)
    except Exception as exc:
        return [*errors, f"Cannot parse metadata TSV: {exc}"]
    if "Protein_ID" not in metadata:
        errors.append("Metadata is missing Protein_ID")
        return errors
    required_metadata = set(METADATA_COLUMNS) - {"cluster_representative_in_reference"}
    missing_columns = required_metadata - set(metadata.columns)
    if missing_columns:
        errors.append("Metadata is missing required columns: " + ", ".join(sorted(missing_columns)))
    if metadata["Protein_ID"].duplicated().any():
        errors.append("Metadata contains duplicate Protein_ID values")
    metadata_ids = set(metadata["Protein_ID"])
    if set(fasta_ids) != metadata_ids:
        errors.append("Metadata Protein_ID values do not exactly match reference FASTA IDs")
    if "Predicted_v3" not in metadata or not set(metadata["Predicted_v3"]) <= {"active", "silent"}:
        errors.append("Metadata Predicted_v3 values are not restricted to active/silent")
    if "Cluster_c4_representative" in metadata:
        if metadata["Cluster_c4_representative"].isin({"", "NA"}).any():
            errors.append("Metadata contains empty Cluster_c4_representative values")
        representatives = set(metadata["Cluster_c4_representative"]) - {"", "NA"}
        missing_representatives = representatives - metadata_ids
        if "cluster_representative_in_reference" in metadata:
            expected = metadata.loc[
                metadata["cluster_representative_in_reference"].str.lower().isin({"true", "1"}),
                "Cluster_c4_representative",
            ]
            missing_representatives = set(expected) - metadata_ids
        else:
            # Older/minimal databases cannot distinguish intentionally excluded
            # unphenotyped representatives. Existing representatives are still checked.
            missing_representatives = set()
        if missing_representatives:
            errors.append(
                "Cluster representatives marked as present but absent from metadata: "
                + ", ".join(sorted(missing_representatives)[:10])
            )
    mmseqs_prefix = database / "flapro_mmseqs" / "refDB"
    mmseqs_files = [mmseqs_prefix, mmseqs_prefix.with_suffix(".dbtype"), mmseqs_prefix.with_suffix(".index")]
    if not all(path.is_file() for path in mmseqs_files):
        errors.append(f"MMseqs2 database files missing for prefix {mmseqs_prefix}")
    return errors
