from __future__ import annotations

import pandas as pd
import pytest

from flapro_genome.classify import Thresholds, classify_hits


def metadata() -> pd.DataFrame:
    rows = [
        ("REF_A", "CL_A", "active", False),
        ("REF_B", "CL_B", "silent", True),
        ("REF_A2", "CL_A", "silent", False),
    ]
    return pd.DataFrame(
        [
            {
                "Protein_ID": protein,
                "Cluster_c4_representative": cluster,
                "Cluster_CDHit": "10",
                "Phylum": "Firmicutes",
                "Class": "Bacilli",
                "Order": "Order",
                "Family": "Family",
                "Genus": "Genus",
                "Species": "Species name",
                "Predicted_v3": phenotype,
                "Experimental": "not_checked",
                "MarkerType": "True Marker",
                "cluster_phenotype_consistent": consistent,
            }
            for protein, cluster, phenotype, consistent in rows
        ]
    )


def hit(
    query: str,
    target: str,
    pident: float,
    qcov: float = 1.0,
    tcov: float = 1.0,
    bitscore: float = 500.0,
    evalue: float = 1e-100,
) -> dict[str, float | str | int]:
    qlen = tlen = 500
    alnlen = int(qcov * qlen)
    return {
        "query": query, "target": target, "pident": pident, "alnlen": alnlen,
        "qlen": qlen, "tlen": tlen, "qcov": qcov, "tcov": tcov,
        "evalue": evalue, "bitscore": bitscore,
    }


@pytest.mark.parametrize(
    ("row", "expected"),
    [
        (hit("Q1", "REF_A", 100), "exact_reference"),
        (hit("Q1", "REF_A", 96), "known_flapro_type"),
        (hit("Q1", "REF_A", 70), "related_flagellin"),
        (hit("Q1", "REF_A", 96, qcov=0.3, tcov=0.6), "weak_candidate"),
    ],
)
def test_classification_levels(row: dict, expected: str) -> None:
    result = classify_hits(pd.DataFrame([row]), metadata(), Thresholds(), "genome")
    assert result.iloc[0]["classification"] == expected


def test_near_best_hits_from_two_clusters_are_ambiguous() -> None:
    hits = pd.DataFrame(
        [
            hit("Q1", "REF_A", 96, bitscore=500),
            hit("Q1", "REF_B", 96, bitscore=495),
        ]
    )
    result = classify_hits(hits, metadata(), Thresholds(), "genome").iloc[0]
    assert result["classification"] == "ambiguous"
    assert pd.isna(result["assigned_cluster"])
    assert result["n_near_best_clusters"] == 2


def test_mixed_cluster_transfers_nearest_individual_phenotype() -> None:
    hits = pd.DataFrame(
        [
            hit("Q1", "REF_A", 100, bitscore=500),
            hit("Q1", "REF_A2", 99.5, bitscore=495),
        ]
    )
    result = classify_hits(hits, metadata(), Thresholds(), "genome").iloc[0]
    assert result["classification"] == "exact_reference"
    assert result["cluster_phenotype_consistent"] == False  # noqa: E712
    assert result["tlr5_class"] == "active"
    assert result["tlr5_reference_protein"] == "REF_A"

