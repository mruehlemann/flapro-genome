from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from flapro_genome.build_db import join_references, parse_parameter_tsv
from flapro_genome.utils import FlaProGenomeError


def test_parse_malformed_parameter_tsv_recovers_protein_id(tmp_path: Path) -> None:
    source = tmp_path / "parameters.tsv"
    source.write_text(
        "Cluster_CDHit\tsequence\tExp.phenotype\n"
        "PROT_A\t12\tMSTNQK\tNA\n"
        "PROT_B\t13\tMKKLAA\tactive\n"
    )
    result = parse_parameter_tsv(source)
    assert result["Protein_ID"].tolist() == ["PROT_A", "PROT_B"]
    assert result["Cluster_CDHit"].tolist() == ["12", "13"]
    assert result["length"].tolist() == [6, 6]


def test_parse_parameter_tsv_rejects_unexpected_layout(tmp_path: Path) -> None:
    source = tmp_path / "bad.tsv"
    source.write_text("Cluster_CDHit\tsequence\nA\t1\tMKK\textra\n")
    with pytest.raises(FlaProGenomeError, match="Unexpected parameter TSV layout"):
        parse_parameter_tsv(source)


def test_reference_join_reports_unmatched_ids() -> None:
    parameters = pd.DataFrame(
        {
            "Protein_ID": ["A", "B"],
            "Cluster_CDHit": ["1", "2"],
            "sequence": ["MKK", "MAA"],
            "length": [3, 3],
        }
    )
    taxonomy = pd.DataFrame(
        {
            "Flagellin_ID": ["A", "C"],
            "Phylum": ["P", "P"], "Class": ["C", "C"], "Order": ["O", "O"],
            "Family": ["F", "F"], "Genus": ["G", "G"], "Species": ["S", "S"],
            "Predicted_v1": ["non-silent", "non-silent"],
            "Experimental": ["not_checked", "not_checked"],
            "Predicted_v3": ["active", "silent"],
            "Cluster_c4_representative": ["A", "C"],
            "MarkerType": ["True Marker", "True Marker"],
            "Content_of_the_cluster": ["A", "C"],
        }
    )
    joined, stats = join_references(parameters, taxonomy)
    assert joined["Protein_ID"].tolist() == ["A"]
    assert stats == {"n_ids_only_parameter": 1, "n_ids_only_taxonomy": 1, "n_ids_both": 1}

