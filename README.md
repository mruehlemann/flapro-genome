# flapro-genome

`flapro-genome` annotates bacterial flagellin proteins encoded by isolate genomes or MAGs using the full-length reference sequences, sequence clusters, taxonomy, and predicted TLR5 phenotypes distributed with [FlaPro](https://github.com/leylabmpi/FlaPro).

This is not FlaPro's ShortBRED metagenomic quantification workflow. It searches every encoded protein against full-length FlaPro flagellins. Detection of a plausible flagellin is deliberately separate from assignment to a known FlaPro sequence type.

## Installation

```bash
git clone <new repository>
cd flapro-genome

conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

conda env create -f environment.yml
conda activate flapro-genome
```

Verify the installation:

```bash
mmseqs version
prodigal -v
hmmsearch -h | head
flapro-genome --help
```

HMMER is installed for optional future domain validation but is not used by the version 1 pipeline. Prodigal is used only for nucleotide inputs.

## Obtain FlaPro and build the database

The source repository is needed only while building the standalone database.

```bash
git clone https://github.com/leylabmpi/FlaPro.git

flapro-genome build-db \
    --flapro-repo FlaPro \
    --out flapro_db

flapro-genome validate-db --db flapro_db
```

The builder parses FlaPro's slightly malformed model parameter table, joins protein IDs to the taxonomy and phenotype table, writes the full-length reference FASTA and metadata, calculates cluster phenotype consistency, and builds an MMseqs2 database. It reports IDs found only on one side of the join rather than silently discarding them.

The database directory contains:

```text
flapro_db/
├── flapro_reference.faa
├── flapro_reference_metadata.tsv
├── flapro_mmseqs/
├── build_stats.json
└── VERSION
```

`VERSION` records the FlaPro commit, UTC build time, flapro-genome version, and MMseqs2 version.

## Annotate one genome

Protein FASTA:

```bash
flapro-genome annotate \
    genome.faa \
    --db flapro_db \
    --out results/genome1 \
    --threads 16
```

Nucleotide genome or MAG FASTA:

```bash
flapro-genome annotate \
    genome.fna \
    --db flapro_db \
    --out results/genome1 \
    --threads 16
```

GenBank input (`.gb`, `.gbf`, `.gbk`, `.gbff`, or `.genbank`):

```bash
flapro-genome annotate \
    genome.gbk \
    --db flapro_db \
    --out results/genome1 \
    --threads 16
```

For GenBank files, annotated CDS translations are used directly. Protein identifiers are selected in this order: `protein_id`, `locus_tag`, `gene`, then a generated record/CDS identifier. CDS translations missing from the file are derived using `codon_start` and `transl_table` (bacterial genetic code 11 by default). Pseudogene features are skipped. If the file contains no usable CDS proteins, its nucleotide records are passed to Prodigal.

Nucleotide FASTA input is detected from the sequence alphabet, not the filename. It is translated with Prodigal in metagenome mode by default. Use `--prodigal-mode single` for an appropriate complete isolate genome. Gzip-compressed FASTA and GenBank input is supported.

The output prefix above produces:

```text
results/genome1.flagellins.tsv
results/genome1.all_hits.tsv
results/genome1.flagellins.faa
results/genome1.summary.tsv
```

When Prodigal is used, its proteins, genes, and GFF are also retained with `.prodigal.*` suffixes. Original protein FASTA identifiers and GenBank CDS identifiers are preserved.

## Batch mode

```bash
flapro-genome batch \
    --input genomes \
    --db flapro_db \
    --out results \
    --threads 32
```

Files ending in `.faa`, `.fasta`, `.fa`, `.fna`, `.gb`, `.gbf`, `.gbk`, `.gbff`, `.genbank`, or their `.gz` forms are processed sequentially. Each genome gets its own result directory. Combined tables are written to `results/all_flagellins.tsv` and `results/all_genome_summaries.tsv`.

## Classification semantics

Default calls are conservative and all numerical thresholds can be changed in `annotate` and `batch` (see `--help`).

- `exact_reference`: at least 99% identity and at least 95% coverage of both query and reference.
- `known_flapro_type`: at least 95% identity and at least 80% coverage of both sequences, with an unambiguous cluster among near-best hits.
- `related_flagellin`: at least 40% identity, at least 60% coverage of both sequences, and E-value at most 1e-10. The nearest FlaPro cluster is reported but not assigned.
- `ambiguous`: near-best hits (at least 98% of the best bitscore) support multiple FlaPro clusters. No cluster or TLR5 phenotype is assigned.
- `weak_candidate`: passes loose detection (E-value at most 1e-5, bitscore at least 50, and at least 50% query or target coverage) but not a stronger classification.

`active` and `silent` are FlaPro's predicted TLR5 phenotypes. Version 1 transfers this phenotype only for `exact_reference` and `known_flapro_type`, and always records the individual reference protein used. A mixed-phenotype cluster does not override the nearest individual reference phenotype. `related_flagellin`, `ambiguous`, and `weak_candidate` remain unresolved; similarity to a moderately related protein is not a new phenotype prediction.

Presence in a genome means encoded potential. It does not imply transcription, translation, assembly, or biological expression.

By default, the candidate FASTA contains exact, known-type, and related flagellins. `--include-weak` also writes weak candidates; ambiguous proteins remain in the annotation and raw-hit audit tables but are not forced into the flagellin FASTA.

## Outputs

`*.flagellins.tsv` has one row per detected candidate, including the nearest individual reference, both coverage values, both FlaPro cluster identifiers, ambiguity statistics, reference taxonomy, TLR5 transfer status, and cluster phenotype consistency.

`*.all_hits.tsv` retains up to 20 MMseqs2 hits per query protein for auditing, including target cluster, species, and `Predicted_v3`.

`*.summary.tsv` describes the encoded flagellin repertoire: total predicted proteins, counts by classification, unique assigned clusters, and resolved/unresolved TLR5 counts.

Missing values are written as `NA`.

## Leave-one-reference-out benchmark

```bash
flapro-genome benchmark \
    --db flapro_db \
    --out benchmark \
    --threads 16
```

The benchmark excludes each self-hit and evaluates identity thresholds 50–99% crossed with bilateral coverage thresholds 0.5–0.9. It writes:

- `threshold_performance.tsv`: recovery coverage, accuracy, and ambiguity at every threshold combination.
- `reference_predictions.tsv`: per-reference results at the current conservative 95% identity / 80% bilateral-coverage starting point.

`cluster_accuracy` is calculated among classifiable, unambiguous predictions. `ambiguous_fraction` is calculated among references with at least one threshold-passing non-self hit. Benchmarking does not automatically change operational defaults.

## Tests

```bash
pytest
```

Tests cover malformed FlaPro TSV parsing, join accounting, all classification levels, cross-cluster ambiguity, mixed-phenotype clusters, protein/nucleotide FASTA and GenBank handling, and MMseqs2 self-classification.

## Scope and future extension

Version 1 is similarity-based. It does not use ShortBRED markers and does not reproduce FlaPro's Random Forest, AlphaFold, or Foldseek analyses. A future version can add an RF-based phenotype predictor for genuinely novel flagellins without changing the current separation between candidate detection, cluster assignment, and phenotype inference.
