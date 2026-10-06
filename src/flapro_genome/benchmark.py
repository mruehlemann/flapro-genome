from __future__ import annotations

import logging
import tempfile
from pathlib import Path

import pandas as pd

from .annotate import load_metadata, run_mmseqs_search
from .build_db import validate_database
from .utils import FlaProGenomeError


LOG = logging.getLogger("flapro_genome")
IDENTITY_THRESHOLDS = (50, 60, 70, 80, 90, 95, 97, 99)
COVERAGE_THRESHOLDS = (0.5, 0.6, 0.7, 0.8, 0.9)


def _predict_at_threshold(
    query_hits: pd.DataFrame,
    metadata_by_id: pd.DataFrame,
    identity: float,
    coverage: float,
    near_best_fraction: float = 0.98,
) -> tuple[str | None, bool, str | None, float | None, float | None, float | None]:
    eligible = query_hits.loc[
        (query_hits["pident"] >= identity)
        & (query_hits["qcov"] >= coverage)
        & (query_hits["tcov"] >= coverage)
    ].sort_values(["bitscore", "evalue"], ascending=[False, True])
    if eligible.empty:
        return None, False, None, None, None, None
    best = eligible.iloc[0]
    near = eligible.loc[eligible["bitscore"] >= best["bitscore"] * near_best_fraction]
    clusters = sorted(set(near["target"].map(metadata_by_id["Cluster_c4_representative"])))
    ambiguous = len(clusters) != 1
    predicted = None if ambiguous else clusters[0]
    return (
        predicted,
        ambiguous,
        str(best["target"]),
        float(best["pident"]),
        float(best["qcov"]),
        float(best["tcov"]),
    )


def run_benchmark(database_dir: str | Path, output_dir: str | Path, threads: int = 1) -> tuple[Path, Path]:
    database = Path(database_dir).resolve()
    output = Path(output_dir).resolve()
    errors = validate_database(database)
    if errors:
        raise FlaProGenomeError("Database validation failed:\n- " + "\n- ".join(errors))
    output.mkdir(parents=True, exist_ok=True)
    metadata = load_metadata(database)
    metadata_by_id = metadata.set_index("Protein_ID")

    with tempfile.TemporaryDirectory(prefix=".benchmark-work-", dir=output) as work_name:
        hits = run_mmseqs_search(
            database / "flapro_reference.faa",
            database,
            Path(work_name),
            threads,
            max_seqs=21,
        )
    hits = hits.loc[hits["query"] != hits["target"]].copy()
    grouped = {str(query): frame for query, frame in hits.groupby("query", sort=False)}
    total = len(metadata)
    performance_rows = []
    default_predictions = []

    for identity in IDENTITY_THRESHOLDS:
        for coverage in COVERAGE_THRESHOLDS:
            n_classifiable = n_ambiguous = n_unambiguous = n_correct = 0
            for row in metadata.itertuples(index=False):
                query_hits = grouped.get(row.Protein_ID, pd.DataFrame(columns=hits.columns))
                prediction = _predict_at_threshold(
                    query_hits, metadata_by_id, identity, coverage
                )
                predicted_cluster, ambiguous, best_target, pident, qcov, tcov = prediction
                if best_target is not None:
                    n_classifiable += 1
                    if ambiguous:
                        n_ambiguous += 1
                    else:
                        n_unambiguous += 1
                        if predicted_cluster == row.Cluster_c4_representative:
                            n_correct += 1
                if identity == 95 and coverage == 0.8:
                    default_predictions.append(
                        {
                            "Protein_ID": row.Protein_ID,
                            "known_cluster": row.Cluster_c4_representative,
                            "predicted_cluster": predicted_cluster,
                            "status": (
                                "unclassifiable" if best_target is None else "ambiguous" if ambiguous
                                else "correct" if predicted_cluster == row.Cluster_c4_representative
                                else "incorrect"
                            ),
                            "nearest_nonself_reference": best_target,
                            "pident": pident,
                            "qcov": qcov,
                            "tcov": tcov,
                        }
                    )
            performance_rows.append(
                {
                    "minimum_identity": identity,
                    "minimum_qcov_tcov": coverage,
                    "n_classifiable": n_classifiable,
                    "coverage_fraction": n_classifiable / total if total else 0.0,
                    "cluster_accuracy": n_correct / n_unambiguous if n_unambiguous else float("nan"),
                    "ambiguous_fraction": n_ambiguous / n_classifiable if n_classifiable else float("nan"),
                }
            )

    performance_path = output / "threshold_performance.tsv"
    predictions_path = output / "reference_predictions.tsv"
    pd.DataFrame(performance_rows).to_csv(performance_path, sep="\t", index=False, na_rep="NA")
    pd.DataFrame(default_predictions).to_csv(predictions_path, sep="\t", index=False, na_rep="NA")
    LOG.info("Wrote leave-one-reference-out benchmark for %d references to %s", total, output)
    return performance_path, predictions_path
