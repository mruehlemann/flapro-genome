from __future__ import annotations

import shutil
import random
from pathlib import Path

import pandas as pd
import pytest

from flapro_genome.annotate import run_mmseqs_search
from flapro_genome.build_db import validate_database
from flapro_genome.classify import Thresholds, classify_hits
from flapro_genome.fasta import FastaRecord, write_fasta
from flapro_genome.utils import run_command


@pytest.mark.integration
@pytest.mark.skipif(shutil.which("mmseqs") is None, reason="MMseqs2 not installed")
def test_reference_fasta_self_classifies_exactly(tmp_path: Path) -> None:
    database = tmp_path / "db"
    (database / "flapro_mmseqs").mkdir(parents=True)
    fasta = database / "flapro_reference.faa"
    sequence = "".join(random.Random(7).choices("ACDEFGHIKLMNPQRSTVWY", k=300))
    write_fasta([FastaRecord("REF_A", sequence)], fasta)
    metadata = pd.DataFrame(
        [{
            "Protein_ID": "REF_A", "length": len(sequence), "Cluster_CDHit": "1",
            "Cluster_c4_representative": "REF_A", "Phylum": "P", "Class": "C",
            "Order": "O", "Family": "F", "Genus": "G", "Species": "S",
            "Predicted_v1": "non-silent", "Experimental": "not_checked",
            "Predicted_v3": "active", "MarkerType": "True Marker",
            "Content_of_the_cluster": "REF_A", "cluster_n_active": 1,
            "cluster_n_silent": 0, "cluster_n_total": 1,
            "cluster_phenotype_consistent": True,
        }]
    )
    metadata.to_csv(database / "flapro_reference_metadata.tsv", sep="\t", index=False)
    (database / "build_stats.json").write_text("{}\n")
    (database / "VERSION").write_text("test=true\n")
    run_command(["mmseqs", "createdb", fasta, database / "flapro_mmseqs" / "refDB"])
    assert validate_database(database) == []
    hits = run_mmseqs_search(fasta, database, tmp_path / "work", threads=1)
    result = classify_hits(hits, metadata.astype(str), Thresholds(), "self").iloc[0]
    assert result["classification"] == "exact_reference"
    assert result["nearest_flapro_protein"] == "REF_A"
    assert result["pident"] == pytest.approx(100.0)
    assert result["qcov"] == pytest.approx(1.0)
    assert result["tcov"] == pytest.approx(1.0)
