import csv
import sys

from wfhe_experiments.config import BenchmarkSpec, ExperimentConfig
from wfhe_experiments.runner import run_experiment


def test_runner_collects_repeated_records(tmp_path):
    script = tmp_path / "fake_benchmark.py"
    script.write_text("print('operation,time_ms')\nprint('encode,1.5')\n", encoding="utf-8")
    config = ExperimentConfig(
        name="unit",
        output_dir=tmp_path / "results",
        benchmarks=[
            BenchmarkSpec(
                name="fake",
                executable=sys.executable,
                arguments=[str(script)],
                repeats=2,
            )
        ],
    )

    result = run_experiment(config, repo_root=tmp_path)

    with result.records_csv.open(encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 2
    assert [row["run_index"] for row in rows] == ["1", "2"]
    assert all(row["operation"] == "encode" for row in rows)
    assert result.manifest_json.is_file()
