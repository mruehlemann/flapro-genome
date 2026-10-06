from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .annotate import annotate_fasta, batch_annotate
from .benchmark import run_benchmark
from .build_db import build_database, validate_database
from .classify import Thresholds
from .utils import FlaProGenomeError, configure_logging


def _add_logging_options(parser: argparse.ArgumentParser, suppress_defaults: bool = False) -> None:
    group = parser.add_mutually_exclusive_group()
    default = argparse.SUPPRESS if suppress_defaults else False
    group.add_argument(
        "--verbose", action="store_true", default=default,
        help="Show commands and diagnostic details",
    )
    group.add_argument(
        "--quiet", action="store_true", default=default,
        help="Suppress routine progress messages",
    )


def _add_threshold_options(parser: argparse.ArgumentParser) -> None:
    defaults = Thresholds()
    group = parser.add_argument_group("candidate detection thresholds")
    group.add_argument("--candidate-evalue", type=float, default=defaults.candidate_evalue)
    group.add_argument("--candidate-bitscore", type=float, default=defaults.candidate_bitscore)
    group.add_argument(
        "--candidate-coverage", type=float, default=defaults.candidate_coverage,
        help="Minimum of max(query coverage, target coverage)",
    )
    group = parser.add_argument_group("exact-reference thresholds")
    group.add_argument("--exact-identity", type=float, default=defaults.exact_identity)
    group.add_argument("--exact-qcov", type=float, default=defaults.exact_qcov)
    group.add_argument("--exact-tcov", type=float, default=defaults.exact_tcov)
    group = parser.add_argument_group("known-type thresholds")
    group.add_argument("--known-identity", type=float, default=defaults.known_identity)
    group.add_argument("--known-qcov", type=float, default=defaults.known_qcov)
    group.add_argument("--known-tcov", type=float, default=defaults.known_tcov)
    group = parser.add_argument_group("related-flagellin thresholds")
    group.add_argument("--related-identity", type=float, default=defaults.related_identity)
    group.add_argument("--related-qcov", type=float, default=defaults.related_qcov)
    group.add_argument("--related-tcov", type=float, default=defaults.related_tcov)
    group.add_argument("--related-evalue", type=float, default=defaults.related_evalue)
    group.add_argument(
        "--near-best-fraction", type=float, default=defaults.near_best_fraction,
        help="Hits with this fraction of the best bitscore are checked for cluster ambiguity",
    )


def _add_annotation_runtime_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--db", required=True, type=Path, help="Built flapro-genome database directory")
    parser.add_argument("--out", required=True, type=Path, help="Output prefix (or directory for batch)")
    parser.add_argument("--threads", type=int, default=1, help="Threads passed to MMseqs2")
    parser.add_argument(
        "--prodigal-mode", choices=("meta", "single"), default="meta",
        help="Prodigal mode for nucleotide input (default: meta)",
    )
    parser.add_argument(
        "--include-weak", action="store_true",
        help="Include weak candidates in the candidate protein FASTA",
    )
    _add_threshold_options(parser)


def _thresholds(args: argparse.Namespace) -> Thresholds:
    return Thresholds(
        candidate_evalue=args.candidate_evalue,
        candidate_bitscore=args.candidate_bitscore,
        candidate_coverage=args.candidate_coverage,
        exact_identity=args.exact_identity,
        exact_qcov=args.exact_qcov,
        exact_tcov=args.exact_tcov,
        known_identity=args.known_identity,
        known_qcov=args.known_qcov,
        known_tcov=args.known_tcov,
        related_identity=args.related_identity,
        related_qcov=args.related_qcov,
        related_tcov=args.related_tcov,
        related_evalue=args.related_evalue,
        near_best_fraction=args.near_best_fraction,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="flapro-genome",
        description="Annotate genome-encoded flagellins using FlaPro full-length references.",
    )
    _add_logging_options(parser)
    subparsers = parser.add_subparsers(dest="command", required=True)

    build = subparsers.add_parser("build-db", help="Build a reference database from a FlaPro clone")
    _add_logging_options(build, suppress_defaults=True)
    build.add_argument("--flapro-repo", required=True, type=Path, help="Path to the FlaPro repository")
    build.add_argument("--out", required=True, type=Path, help="New database directory")

    validate = subparsers.add_parser("validate-db", help="Check database integrity")
    _add_logging_options(validate, suppress_defaults=True)
    validate.add_argument("--db", required=True, type=Path, help="Built database directory")

    annotate = subparsers.add_parser("annotate", help="Annotate one protein or nucleotide FASTA")
    _add_logging_options(annotate, suppress_defaults=True)
    annotate.add_argument("input", type=Path, help="Protein or nucleotide FASTA, optionally gzip-compressed")
    _add_annotation_runtime_options(annotate)

    batch = subparsers.add_parser("batch", help="Annotate supported FASTA files in a directory")
    _add_logging_options(batch, suppress_defaults=True)
    batch.add_argument("--input", required=True, type=Path, help="Directory containing FASTA files")
    _add_annotation_runtime_options(batch)

    benchmark = subparsers.add_parser(
        "benchmark", help="Run leave-one-reference-out cluster-recovery benchmarking"
    )
    _add_logging_options(benchmark, suppress_defaults=True)
    benchmark.add_argument("--db", required=True, type=Path, help="Built database directory")
    benchmark.add_argument("--out", required=True, type=Path, help="Benchmark output directory")
    benchmark.add_argument("--threads", type=int, default=1, help="Threads passed to MMseqs2")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    configure_logging(args.verbose, args.quiet)
    try:
        if args.command == "build-db":
            stats = build_database(args.flapro_repo, args.out)
            print(json.dumps(stats, indent=2))
        elif args.command == "validate-db":
            errors = validate_database(args.db)
            if errors:
                for error in errors:
                    print(f"ERROR: {error}", file=sys.stderr)
                return 1
            print(f"Database is valid: {args.db}")
        elif args.command == "annotate":
            annotate_fasta(
                args.input, args.db, args.out, _thresholds(args), args.threads,
                args.prodigal_mode, args.include_weak,
            )
        elif args.command == "batch":
            batch_annotate(
                args.input, args.db, args.out, _thresholds(args), args.threads,
                args.prodigal_mode, args.include_weak,
            )
        elif args.command == "benchmark":
            run_benchmark(args.db, args.out, args.threads)
        else:  # pragma: no cover
            parser.error(f"Unknown command: {args.command}")
    except (FlaProGenomeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
