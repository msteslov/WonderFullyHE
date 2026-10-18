import csv

import pandas as pd

from wfhe_experiments.analysis import generate_report, load_records, observations, summarize


def test_analysis_creates_summary_report_and_plot(tmp_path):
    records = tmp_path / "records.csv"
    with records.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=[
                "benchmark",
                "profile",
                "operation",
                "time_ms",
                "max_abs_error",
                "mean_abs_error",
            ],
        )
        writer.writeheader()
        writer.writerows(
            [
                {
                    "benchmark": "ckks",
                    "profile": "fast",
                    "operation": "encode",
                    "time_ms": 1,
                    "max_abs_error": 1e-6,
                    "mean_abs_error": 1e-7,
                },
                {
                    "benchmark": "ckks",
                    "profile": "fast",
                    "operation": "encode",
                    "time_ms": 3,
                    "max_abs_error": 3e-6,
                    "mean_abs_error": 3e-7,
                },
            ]
        )

    frame = load_records(records)
    summary = summarize(frame)
    report, summary_file, plots = generate_report(records, tmp_path / "report")

    assert summary.loc[summary["metric"] == "time_ms", "mean"].iloc[0] == 2
    assert report.is_file()
    assert summary_file.is_file()
    assert [path.name for path in plots] == ["ckks_latency.png", "ckks_error.png"]


def test_observations_describe_warm_regression_as_slower():
    frame = pd.DataFrame(
        [
            {
                "benchmark": "cts",
                "record_type": "phase",
                "measurement": "apply_first",
                "value": 100.0,
            },
            {
                "benchmark": "cts",
                "record_type": "phase",
                "measurement": "apply_warm",
                "value": 110.0,
            },
        ]
    )

    result = observations(frame)

    assert any("10.0% slower" in item for item in result)
