"""Subprocess runner with provenance and machine-readable outputs."""

import csv
import json
import platform
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

from .config import ExperimentConfig
from .parsers import parse_coeff_to_slot_output, parse_csv_output


@dataclass(frozen=True)
class RunResult:
    """Files produced by one complete experiment."""

    directory: Path
    records_csv: Path
    records_json: Path
    manifest_json: Path


class BenchmarkFailed(RuntimeError):
    """Raised when an invocation exits unsuccessfully without permission."""


def _git_value(repo_root: Path, arguments: list[str]) -> Optional[str]:
    process = subprocess.run(
        ["git", *arguments],
        cwd=repo_root,
        capture_output=True,
        text=True,
        check=False,
    )
    return process.stdout.strip() if process.returncode == 0 else None


def _resolve_executable(repo_root: Path, build_dir: Path, executable: str) -> Path:
    path = Path(executable)
    if not path.is_absolute():
        path = repo_root / build_dir / path
    path = path.resolve()
    if not path.is_file():
        raise FileNotFoundError(f"benchmark executable not found: {path}")
    return path


def _write_records(records: list[dict[str, Any]], csv_path: Path, json_path: Path) -> None:
    if not records:
        raise ValueError("experiment produced no records")

    metadata = ["benchmark", "run_index", "return_code", "wall_time_ms"]
    fields = metadata + sorted({key for row in records for key in row if key not in metadata})
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="raise")
        writer.writeheader()
        writer.writerows(records)
    with json_path.open("w", encoding="utf-8") as stream:
        json.dump(records, stream, ensure_ascii=False, indent=2)
        stream.write("\n")


def run_experiment(
    config: ExperimentConfig,
    repo_root: Path,
    output_dir: Optional[Path] = None,
) -> RunResult:
    """Run all configured benchmarks and persist raw output plus provenance."""

    repo_root = repo_root.resolve()
    started = datetime.now(timezone.utc)
    base_output = output_dir or config.output_dir
    if not base_output.is_absolute():
        base_output = repo_root / base_output
    result_dir = base_output / f"{started.strftime('%Y%m%dT%H%M%SZ')}-{config.name}"
    raw_dir = result_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=False)

    records: list[dict[str, Any]] = []
    invocations: list[dict[str, Any]] = []
    pending_failure: Optional[str] = None

    for benchmark in config.benchmarks:
        executable = _resolve_executable(repo_root, config.build_dir, benchmark.executable)
        for run_index in range(1, benchmark.repeats + 1):
            command = [str(executable), *benchmark.arguments]
            before = time.perf_counter()
            process = subprocess.run(
                command,
                cwd=repo_root,
                capture_output=True,
                text=True,
                timeout=benchmark.timeout_seconds,
                check=False,
            )
            wall_time_ms = (time.perf_counter() - before) * 1000.0
            stem = f"{benchmark.name}-{run_index:03d}"
            (raw_dir / f"{stem}.stdout.csv").write_text(process.stdout, encoding="utf-8")
            (raw_dir / f"{stem}.stderr.txt").write_text(process.stderr, encoding="utf-8")

            invocation = {
                "benchmark": benchmark.name,
                "run_index": run_index,
                "command": command,
                "return_code": process.returncode,
                "wall_time_ms": round(wall_time_ms, 3),
            }
            invocations.append(invocation)

            parsed = (
                parse_coeff_to_slot_output(process.stdout)
                if benchmark.parser == "coeff_to_slot"
                else parse_csv_output(process.stdout)
            )
            for row in parsed:
                records.append(
                    {
                        "benchmark": benchmark.name,
                        "run_index": run_index,
                        "return_code": process.returncode,
                        "wall_time_ms": round(wall_time_ms, 3),
                        **row,
                    }
                )

            if process.returncode != 0 and not benchmark.allow_nonzero:
                pending_failure = (
                    f"{benchmark.name} run {run_index} exited with code {process.returncode}"
                )
                break
        if pending_failure:
            break

    records_csv = result_dir / "records.csv"
    records_json = result_dir / "records.json"
    manifest_json = result_dir / "manifest.json"
    _write_records(records, records_csv, records_json)

    status = _git_value(repo_root, ["status", "--porcelain"])
    manifest = {
        "schema_version": 1,
        "experiment": config.model_dump(mode="json"),
        "started_at_utc": started.isoformat(),
        "finished_at_utc": datetime.now(timezone.utc).isoformat(),
        "repository": {
            "root": str(repo_root),
            "commit": _git_value(repo_root, ["rev-parse", "HEAD"]),
            "dirty": bool(status) if status is not None else None,
        },
        "runtime": {
            "python": sys.version,
            "platform": platform.platform(),
        },
        "invocations": invocations,
    }
    with manifest_json.open("w", encoding="utf-8") as stream:
        json.dump(manifest, stream, ensure_ascii=False, indent=2)
        stream.write("\n")

    result = RunResult(result_dir, records_csv, records_json, manifest_json)
    if pending_failure:
        raise BenchmarkFailed(f"{pending_failure}; partial results saved to {result_dir}")
    return result
