"""Command-line interface for reproducible WonderFullyHE experiments."""

import argparse
from collections.abc import Sequence
from pathlib import Path
from typing import Optional

from .analysis import generate_report
from .config import load_config
from .runner import run_experiment


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="wfhe-experiments")
    subparsers = parser.add_subparsers(dest="command", required=True)

    run = subparsers.add_parser("run", help="run benchmarks from a JSON config")
    run.add_argument("config", type=Path)
    run.add_argument("--repo-root", type=Path, default=Path.cwd())
    run.add_argument("--output-dir", type=Path)
    run.add_argument("--no-report", action="store_true")

    report = subparsers.add_parser("report", help="analyze an existing records.csv")
    report.add_argument("records", type=Path)
    report.add_argument("--output-dir", type=Path)
    return parser


def main(arguments: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(arguments)
    if args.command == "run":
        result = run_experiment(
            load_config(args.config),
            repo_root=args.repo_root,
            output_dir=args.output_dir,
        )
        if not args.no_report:
            generate_report(result.records_csv, result.directory)
        print(result.directory)
        return 0

    output_dir = args.output_dir or args.records.parent
    report_path, _, _ = generate_report(args.records, output_dir)
    print(report_path)
    return 0
