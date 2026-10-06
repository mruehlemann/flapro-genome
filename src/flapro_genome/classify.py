from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import pandas as pd


@dataclass(frozen=True)
class Thresholds:
    candidate_evalue: float = 1e-5
    candidate_bitscore: float = 50.0
    candidate_coverage: float = 0.50
    exact_identity: float = 99.0
    exact_qcov: float = 0.95
    exact_tcov: float = 0.95
    known_identity: float = 95.0
    known_qcov: float = 0.80
    known_tcov: float = 0.80
    related_identity: float = 40.0
    related_qcov: float = 0.60
    related_tcov: float = 0.60
    related_evalue: float = 1e-10
    near_best_fraction: float = 0.98

    def __post_init__(self) -> None:
        for name in ("candidate_evalue", "related_evalue"):
            if getattr(self, name) <= 0:
                raise ValueError(f"{name} must be greater than zero")
        if self.candidate_bitscore < 0:
            raise ValueError("candidate_bitscore must be non-negative")
        for name in (
            "candidate_coverage", "exact_qcov", "exact_tcov", "known_qcov", "known_tcov",
            "related_qcov", "related_tcov", "near_best_fraction",
        ):
            value = getattr(self, name)
            if not 0 <= value <= 1:
                raise ValueError(f"{name} must be between 0 and 1")
        for name in ("exact_identity", "known_identity", "related_identity"):
            value = getattr(self, name)
            if not 0 <= value <= 100:
                raise ValueError(f"{name} must be between 0 and 100")


CLASSIFICATION_ORDER = {
    "exact_reference": 0,
    "known_flapro_type": 1,
    "related_flagellin": 2,
    "ambiguous": 3,
    "weak_candidate": 4,
}


def candidate_mask(hits: pd.DataFrame, thresholds: Thresholds) -> pd.Series:
    return (
        (hits["evalue"] <= thresholds.candidate_evalue)
        & (hits["bitscore"] >= thresholds.candidate_bitscore)
        & (hits[["qcov", "tcov"]].max(axis=1) >= thresholds.candidate_coverage)
    )


def _passes(best: pd.Series, identity: float, qcov: float, tcov: float, evalue: float | None = None) -> bool:
    passed = best["pident"] >= identity and best["qcov"] >= qcov and best["tcov"] >= tcov
    return bool(passed and (evalue is None or best["evalue"] <= evalue))


def _value(metadata: pd.Series, key: str) -> Any:
    value = metadata.get(key, pd.NA)
    if pd.isna(value) or value == "" or value == "NA":
        return pd.NA
    return value


def classify_query_hits(
    query: str,
    hits: pd.DataFrame,
    metadata_by_id: pd.DataFrame,
    thresholds: Thresholds,
    genome: str,
) -> dict[str, Any] | None:
    """Classify one query using candidate-qualified hits and near-best ambiguity."""
    qualifying = hits.loc[candidate_mask(hits, thresholds)].copy()
    qualifying["_is_self_id"] = qualifying["target"].astype(str).eq(str(query))
    qualifying = qualifying.sort_values(
        ["bitscore", "evalue", "pident", "_is_self_id"],
        ascending=[False, True, False, False],
    )
    if qualifying.empty:
        return None
    best = qualifying.iloc[0]
    target = str(best["target"])
    target_metadata = metadata_by_id.loc[target]
    if isinstance(target_metadata, pd.DataFrame):
        raise ValueError(f"Reference metadata contains duplicate ID {target}")

    near_best = qualifying.loc[qualifying["bitscore"] >= float(best["bitscore"]) * thresholds.near_best_fraction]
    clusters = sorted(
        {
            str(metadata_by_id.loc[str(target_id), "Cluster_c4_representative"])
            for target_id in near_best["target"]
            if str(metadata_by_id.loc[str(target_id), "Cluster_c4_representative"]) not in {"", "NA", "nan"}
        }
    )
    cluster_consistent = len(clusters) <= 1
    nearest_cluster = _value(target_metadata, "Cluster_c4_representative")

    if _passes(best, thresholds.exact_identity, thresholds.exact_qcov, thresholds.exact_tcov):
        preliminary = "exact_reference"
    elif _passes(best, thresholds.known_identity, thresholds.known_qcov, thresholds.known_tcov):
        preliminary = "known_flapro_type"
    elif _passes(
        best,
        thresholds.related_identity,
        thresholds.related_qcov,
        thresholds.related_tcov,
        thresholds.related_evalue,
    ):
        preliminary = "related_flagellin"
    else:
        preliminary = "weak_candidate"

    classification = preliminary
    if not cluster_consistent:
        classification = "ambiguous"
    assigned_cluster: Any = pd.NA
    if classification in {"exact_reference", "known_flapro_type"} and pd.notna(nearest_cluster):
        assigned_cluster = nearest_cluster

    transfer = classification in {"exact_reference", "known_flapro_type"}
    nearest_tlr5 = _value(target_metadata, "Predicted_v3")
    row: dict[str, Any] = {
        "genome": genome,
        "query_protein": query,
        "query_length": int(best["qlen"]),
        "classification": classification,
        "nearest_flapro_protein": target,
        "pident": float(best["pident"]),
        "alignment_length": int(best["alnlen"]),
        "qcov": float(best["qcov"]),
        "tcov": float(best["tcov"]),
        "evalue": float(best["evalue"]),
        "bitscore": float(best["bitscore"]),
        "assigned_cluster": assigned_cluster,
        "nearest_cluster": nearest_cluster,
        "cluster_consistent": cluster_consistent,
        "n_near_best_hits": int(len(near_best)),
        "n_near_best_clusters": int(len(clusters)),
        "near_best_clusters": ",".join(clusters) if clusters else pd.NA,
        "Cluster_CDHit": _value(target_metadata, "Cluster_CDHit"),
        "reference_phylum": _value(target_metadata, "Phylum"),
        "reference_class": _value(target_metadata, "Class"),
        "reference_order": _value(target_metadata, "Order"),
        "reference_family": _value(target_metadata, "Family"),
        "reference_genus": _value(target_metadata, "Genus"),
        "reference_species": _value(target_metadata, "Species"),
        "tlr5_class": nearest_tlr5 if transfer else pd.NA,
        "tlr5_source": "reference_transfer" if transfer else "none",
        "nearest_reference_tlr5_class": nearest_tlr5,
        "cluster_phenotype_consistent": _value(target_metadata, "cluster_phenotype_consistent"),
        "tlr5_reference_protein": target if transfer else pd.NA,
        "tlr5_reference_identity": float(best["pident"]) if transfer else pd.NA,
        "reference_experimental_phenotype": _value(target_metadata, "Experimental"),
        "reference_marker_type": _value(target_metadata, "MarkerType"),
    }
    return row


def classify_hits(
    hits: pd.DataFrame,
    metadata: pd.DataFrame,
    thresholds: Thresholds,
    genome: str,
) -> pd.DataFrame:
    metadata_by_id = metadata.set_index("Protein_ID", drop=False)
    rows = []
    for query, query_hits in hits.groupby("query", sort=False):
        row = classify_query_hits(str(query), query_hits, metadata_by_id, thresholds, genome)
        if row is not None:
            rows.append(row)
    result = pd.DataFrame(rows)
    if not result.empty:
        result["_rank"] = result["classification"].map(CLASSIFICATION_ORDER)
        result = result.sort_values(["_rank", "bitscore"], ascending=[True, False]).drop(columns="_rank")
    return result.reset_index(drop=True)
