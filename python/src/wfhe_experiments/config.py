"""Validated configuration models for benchmark experiments."""

import json
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, Field, field_validator


class BenchmarkSpec(BaseModel):
    """One executable invocation repeated a fixed number of times."""

    name: str = Field(min_length=1, pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
    executable: str = Field(min_length=1)
    arguments: list[str] = Field(default_factory=list)
    repeats: int = Field(default=1, ge=1, le=1000)
    timeout_seconds: float = Field(default=300.0, gt=0)
    allow_nonzero: bool = False
    parser: Literal["csv", "coeff_to_slot"] = "csv"


class ExperimentConfig(BaseModel):
    """Top-level experiment description."""

    name: str = Field(min_length=1, pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
    build_dir: Path = Path("build/bin")
    output_dir: Path = Path("python/results")
    benchmarks: list[BenchmarkSpec] = Field(min_length=1)

    @field_validator("benchmarks")
    @classmethod
    def benchmark_names_must_be_unique(
        cls, benchmarks: list[BenchmarkSpec]
    ) -> list[BenchmarkSpec]:
        names = [item.name for item in benchmarks]
        if len(names) != len(set(names)):
            raise ValueError("benchmark names must be unique")
        return benchmarks


def load_config(path: Path) -> ExperimentConfig:
    """Load and validate a JSON experiment configuration."""

    with path.open(encoding="utf-8") as stream:
        return ExperimentConfig.model_validate(json.load(stream))
