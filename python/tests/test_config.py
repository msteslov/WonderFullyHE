import json

import pytest
from pydantic import ValidationError

from wfhe_experiments.config import load_config


def test_load_config(tmp_path):
    path = tmp_path / "experiment.json"
    path.write_text(
        json.dumps(
            {
                "name": "test",
                "benchmarks": [{"name": "bench", "executable": "binary", "repeats": 2}],
            }
        ),
        encoding="utf-8",
    )

    config = load_config(path)

    assert config.name == "test"
    assert config.benchmarks[0].repeats == 2


def test_duplicate_benchmark_names_are_rejected(tmp_path):
    path = tmp_path / "experiment.json"
    path.write_text(
        json.dumps(
            {
                "name": "test",
                "benchmarks": [
                    {"name": "same", "executable": "first"},
                    {"name": "same", "executable": "second"},
                ],
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(ValidationError, match="benchmark names must be unique"):
        load_config(path)
