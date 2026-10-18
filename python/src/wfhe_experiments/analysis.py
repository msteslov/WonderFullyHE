"""Statistical summaries and plots for experiment records."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.switch_backend("Agg")


NUMERIC_COLUMNS = {
    "run_index",
    "return_code",
    "wall_time_ms",
    "polyModulusDegree",
    "payload_size",
    "time_ms",
    "max_abs_error",
    "mean_abs_error",
    "serialized_bytes",
    "active_slots",
    "poly_modulus_degree",
    "scale_log2",
    "trials",
    "key_setup_ms",
    "add_max_abs_error",
    "multiply_max_abs_error",
    "rotate_max_abs_error",
    "max_rotate_ms",
    "input_chain_index",
    "multiply_chain_index",
    "public_key_bytes",
    "relin_key_bytes",
    "galois_key_bytes",
    "value",
    "work_bits",
    "chain_length",
    "work_levels",
    "depth",
    "mean_max_abs_error",
    "worst_max_abs_error",
    "mean_mean_abs_error",
    "worst_mean_abs_error",
    "threads",
    "tasks",
    "setup_ms",
    "runtime_ms",
    "throughput_tasks_per_sec",
    "tolerance",
    "estimated_diagonals",
    "estimated_rotations",
    "estimated_keys",
    "actual_multiplications",
    "actual_rotations",
    "actual_keys",
    "context_ms",
    "keys_ms",
    "input_ms",
    "plan_ms",
    "apply_ms",
    "decode_ms",
    "total_ms",
    "nodes",
    "levels_used",
    "output_chain_index",
    "observed_error",
    "target_error",
    "test_only",
}

SUMMARY_METRICS = (
    "time_ms",
    "wall_time_ms",
    "key_setup_ms",
    "max_rotate_ms",
    "mean_max_abs_error",
    "worst_max_abs_error",
    "mean_mean_abs_error",
    "worst_mean_abs_error",
    "setup_ms",
    "runtime_ms",
    "throughput_tasks_per_sec",
    "max_abs_error",
    "mean_abs_error",
    "value",
    "context_ms",
    "keys_ms",
    "input_ms",
    "plan_ms",
    "apply_ms",
    "decode_ms",
    "total_ms",
    "observed_error",
)

DIMENSIONS = (
    "benchmark",
    "profile",
    "operation",
    "candidate",
    "experiment",
    "variant",
    "threads",
    "record_type",
    "measurement",
    "unit",
    "scenario",
    "test_only",
)


def load_records(path: Path) -> pd.DataFrame:
    """Load records and coerce known measurement columns to numeric values."""

    frame = pd.read_csv(path)
    for column in NUMERIC_COLUMNS.intersection(frame.columns):
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    return frame


def summarize(frame: pd.DataFrame) -> pd.DataFrame:
    """Calculate long-form statistics and normal-approximation 95% intervals."""

    summaries: list[pd.DataFrame] = []
    for metric in SUMMARY_METRICS:
        if metric not in frame.columns:
            continue
        rows = frame.loc[frame[metric].notna()].copy()
        if rows.empty:
            continue
        dimensions = [
            column for column in DIMENSIONS if column in rows.columns and rows[column].notna().any()
        ]
        if not dimensions:
            continue
        statistics = rows.groupby(dimensions, dropna=False)[metric].agg(
            ["count", "mean", "median", "min", "max", "std"]
        )
        statistics["ci95"] = 1.96 * statistics["std"] / np.sqrt(statistics["count"])
        statistics.insert(0, "metric", metric)
        summaries.append(statistics.reset_index())

    if not summaries:
        return pd.DataFrame()
    return pd.concat(summaries, ignore_index=True, sort=False)


def _markdown_table(frame: pd.DataFrame, maximum_rows: int = 80) -> str:
    if frame.empty:
        return "Нет данных для сводки."
    shown = frame.head(maximum_rows).fillna("")
    columns = list(shown.columns)
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for _, row in shown.iterrows():
        values = [str(row[column]).replace("|", "\\|") for column in columns]
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def _save_ckks_plot(frame: pd.DataFrame, output_dir: Path) -> list[Path]:
    required = {"operation", "profile", "time_ms"}
    subset = (
        frame.dropna(subset=list(required))
        if required.issubset(frame.columns)
        else frame.iloc[0:0]
    )
    if subset.empty:
        return []

    latency = subset.groupby(["operation", "profile"])["time_ms"].mean().unstack()
    axis = latency.plot(kind="bar", figsize=(12, 6))
    axis.set_title("CKKS operation latency")
    axis.set_ylabel("Mean time, ms")
    axis.set_xlabel("Operation")
    axis.tick_params(axis="x", rotation=35)
    axis.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    path = output_dir / "ckks_latency.png"
    plt.savefig(path, dpi=160)
    plt.close()

    plots = [path]
    error_columns = [
        column
        for column in ("max_abs_error", "mean_abs_error")
        if column in subset.columns and subset[column].notna().any()
    ]
    if error_columns:
        errors = subset.groupby("operation")[error_columns].mean()
        errors = errors[(errors > 0).any(axis=1)]
        if not errors.empty:
            axis = errors.plot(kind="bar", figsize=(12, 6), logy=True)
            axis.set_title("CKKS numerical error")
            axis.set_ylabel("Mean absolute error, logarithmic scale")
            axis.set_xlabel("Operation")
            axis.tick_params(axis="x", rotation=35)
            axis.grid(axis="y", alpha=0.25)
            plt.tight_layout()
            error_path = output_dir / "ckks_error.png"
            plt.savefig(error_path, dpi=160)
            plt.close()
            plots.append(error_path)
    return plots


def _save_bootstrap_plot(frame: pd.DataFrame, output_dir: Path) -> list[Path]:
    required = {"candidate", "key_setup_ms"}
    subset = (
        frame.dropna(subset=list(required))
        if required.issubset(frame.columns)
        else frame.iloc[0:0]
    )
    if subset.empty:
        return []

    values = subset.groupby("candidate")["key_setup_ms"].mean().sort_values()
    axis = values.plot(kind="bar", figsize=(10, 5), color="#476f95")
    axis.set_title("Bootstrap candidate key setup")
    axis.set_ylabel("Mean time, ms")
    axis.set_xlabel("Candidate")
    axis.tick_params(axis="x", rotation=25)
    axis.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    path = output_dir / "bootstrap_key_setup.png"
    plt.savefig(path, dpi=160)
    plt.close()
    return [path]


def _save_chain_accuracy_plot(frame: pd.DataFrame, output_dir: Path) -> list[Path]:
    required = {"experiment", "variant", "work_levels", "worst_max_abs_error"}
    if not required.issubset(frame.columns):
        return []
    subset = frame.loc[
        (frame["experiment"] == "chain_length_fixed_depth")
        & (frame.get("status", "ok") == "ok")
    ].dropna(subset=list(required))
    if subset.empty:
        return []

    values = (
        subset.groupby(["work_levels", "variant"])["worst_max_abs_error"].mean().unstack()
    )
    axis = values.plot(marker="o", figsize=(10, 5), logy=True)
    axis.set_title("CKKS error by multiplicative depth and chain length")
    axis.set_ylabel("Worst maximum absolute error")
    axis.set_xlabel("Work levels")
    axis.grid(alpha=0.25)
    plt.tight_layout()
    path = output_dir / "chain_accuracy.png"
    plt.savefig(path, dpi=160)
    plt.close()
    return [path]


def _save_throughput_plot(frame: pd.DataFrame, output_dir: Path) -> list[Path]:
    required = {"threads", "throughput_tasks_per_sec"}
    if not required.issubset(frame.columns):
        return []
    subset = frame.dropna(subset=list(required))
    if "status" in subset.columns:
        subset = subset.loc[subset["status"].str.upper().isin({"PASS", "OK"})]
    if subset.empty:
        return []

    values = subset.groupby("threads")["throughput_tasks_per_sec"].mean().sort_index()
    axis = values.plot(marker="o", figsize=(8, 5), color="#287271")
    axis.set_title("Parallel CKKS throughput")
    axis.set_ylabel("Tasks per second")
    axis.set_xlabel("Worker threads")
    axis.grid(alpha=0.25)
    plt.tight_layout()
    path = output_dir / "parallel_throughput.png"
    plt.savefig(path, dpi=160)
    plt.close()
    return [path]


def _save_coeff_to_slot_plot(frame: pd.DataFrame, output_dir: Path) -> list[Path]:
    required = {"record_type", "measurement", "value"}
    if not required.issubset(frame.columns):
        return []
    subset = frame.loc[frame["record_type"] == "phase"].dropna(subset=list(required))
    if subset.empty:
        return []

    values = subset.groupby("measurement")["value"].mean().sort_values()
    axis = values.plot(kind="barh", figsize=(9, 6), color="#8f5d5d")
    axis.set_title("CoeffToSlot phase latency")
    axis.set_xlabel("Mean time, ms")
    axis.set_ylabel("Phase")
    axis.grid(axis="x", alpha=0.25)
    plt.tight_layout()
    path = output_dir / "coeff_to_slot_phases.png"
    plt.savefig(path, dpi=160)
    plt.close()
    return [path]


def _save_end_to_end_plot(frame: pd.DataFrame, output_dir: Path) -> list[Path]:
    metrics = [
        column
        for column in ("context_ms", "keys_ms", "input_ms", "plan_ms", "apply_ms", "decode_ms")
        if column in frame.columns and frame[column].notna().any()
    ]
    if "scenario" not in frame.columns or not metrics:
        return []
    subset = frame.loc[frame["scenario"] == "test_fixture"]
    if subset.empty:
        return []

    values = subset[metrics].mean().sort_values()
    axis = values.plot(kind="barh", figsize=(9, 5), color="#665191")
    axis.set_title("Test-only end-to-end bootstrap latency")
    axis.set_xlabel("Mean time, ms")
    axis.set_ylabel("Phase")
    axis.grid(axis="x", alpha=0.25)
    plt.tight_layout()
    path = output_dir / "bootstrap_end_to_end_fixture.png"
    plt.savefig(path, dpi=160)
    plt.close()
    return [path]


def observations(frame: pd.DataFrame) -> list[str]:
    """Build cautious, data-derived observations without claiming causality."""

    result: list[str] = []
    if {"profile", "operation", "time_ms", "max_abs_error"}.issubset(frame.columns):
        multiplication = frame.loc[
            frame["operation"] == "multiplyRelinearizeAndRescale"
        ].dropna(subset=["profile", "time_ms", "max_abs_error"])
        profile_means = multiplication.groupby("profile")[["time_ms", "max_abs_error"]].mean()
        if {"basic_ckks", "high_precision_ckks"}.issubset(profile_means.index):
            basic = profile_means.loc["basic_ckks"]
            precision = profile_means.loc["high_precision_ckks"]
            error_improvement = basic["max_abs_error"] / precision["max_abs_error"]
            latency_cost = precision["time_ms"] / basic["time_ms"]
            result.append(
                f"For multiply-rescale, high_precision_ckks reduced mean maximum error "
                f"by {error_improvement:.0f}x versus basic_ckks at {latency_cost:.2f}x latency."
            )

    if {"experiment", "variant", "work_levels", "status"}.issubset(frame.columns):
        chain_failures = frame.loc[
            (frame["experiment"] == "chain_length_fixed_depth")
            & (frame["status"] == "failed")
        ]
        for variant, rows in chain_failures.groupby("variant"):
            levels = sorted({int(level) for level in rows["work_levels"].dropna()})
            if levels:
                rendered_levels = ", ".join(str(level) for level in levels)
                result.append(
                    f"{variant} failed with {rendered_levels} work level(s), indicating an "
                    "insufficient modulus chain for that depth."
                )

    if {"threads", "throughput_tasks_per_sec"}.issubset(frame.columns):
        throughput_rows = frame.dropna(subset=["threads", "throughput_tasks_per_sec"])
        if "status" in throughput_rows.columns:
            failed = throughput_rows.loc[
                ~throughput_rows["status"].str.upper().isin({"PASS", "OK"})
            ]
            if not failed.empty:
                result.append(
                    f"Parallel benchmark rejected {len(failed)}/{len(throughput_rows)} rows "
                    "on its accuracy criterion; rejected rows are excluded from "
                    "throughput comparison."
                )
            throughput_rows = throughput_rows.loc[
                throughput_rows["status"].str.upper().isin({"PASS", "OK"})
            ]
        throughput = throughput_rows.groupby("threads")["throughput_tasks_per_sec"].mean()
        if not throughput.empty and 1 in throughput.index:
            best_threads = throughput.idxmax()
            speedup = throughput.max() / throughput.loc[1]
            result.append(
                f"Maximum observed throughput used {best_threads:g} threads: "
                f"{throughput.max():.2f} tasks/s ({speedup:.2f}x versus one thread)."
            )

    if {"record_type", "measurement", "value"}.issubset(frame.columns):
        phases = frame.loc[frame["record_type"] == "phase"].groupby("measurement")["value"].mean()
        if {"apply_first", "apply_warm"}.issubset(phases.index):
            change = 100.0 * (phases["apply_first"] - phases["apply_warm"]) / phases["apply_first"]
            comparison = "faster" if change >= 0.0 else "slower"
            result.append(
                f"Warm CoeffToSlot application was {abs(change):.1f}% {comparison} than "
                "the first application when averaged across tested factorizations."
            )
        warm_rows = frame.loc[
            (frame["record_type"] == "phase") & (frame["measurement"] == "apply_warm")
        ].dropna(subset=["benchmark", "value"])
        warm_candidates = warm_rows.groupby("benchmark")["value"].mean()
        if not warm_candidates.empty:
            result.append(
                f"Fastest observed warm CoeffToSlot factorization: "
                f"{warm_candidates.idxmin()} at {warm_candidates.min():.2f} ms."
            )

    if {"candidate", "key_setup_ms"}.issubset(frame.columns):
        candidates = frame.dropna(subset=["candidate", "key_setup_ms"]).groupby("candidate")[
            "key_setup_ms"
        ].mean()
        if not candidates.empty:
            result.append(
                f"Lowest measured bootstrap key-setup time: {candidates.idxmin()} "
                f"at {candidates.min():.2f} ms."
            )
    if {"scenario", "apply_ms", "observed_error", "status"}.issubset(frame.columns):
        fixture = frame.loc[frame["scenario"] == "test_fixture"].dropna(
            subset=["apply_ms", "observed_error", "status"]
        )
        if not fixture.empty:
            passed = int((fixture["status"] == "PASS").sum())
            result.append(
                f"Test-only end-to-end fixture passed {passed}/{len(fixture)} runs with "
                f"mean apply latency {fixture['apply_ms'].mean():.2f} ms and mean observed "
                f"error {fixture['observed_error'].mean():.3e}; this is not production evidence."
            )
    return result


def generate_report(records_path: Path, output_dir: Path) -> tuple[Path, Path, list[Path]]:
    """Create a CSV summary, Markdown report, and available plots."""

    output_dir.mkdir(parents=True, exist_ok=True)
    frame = load_records(records_path)
    summary = summarize(frame)
    summary_path = output_dir / "summary.csv"
    summary.to_csv(summary_path, index=False)

    plots = (
        _save_ckks_plot(frame, output_dir)
        + _save_bootstrap_plot(frame, output_dir)
        + _save_chain_accuracy_plot(frame, output_dir)
        + _save_throughput_plot(frame, output_dir)
        + _save_coeff_to_slot_plot(frame, output_dir)
        + _save_end_to_end_plot(frame, output_dir)
    )
    report_path = output_dir / "report.md"
    report = [
        "# WonderFullyHE experiment report",
        "",
        f"Source: `{records_path}`",
        "",
        f"Records: **{len(frame)}**",
        "",
        "## Timing summary",
        "",
        _markdown_table(summary),
    ]
    generated_observations = observations(frame)
    if generated_observations:
        report.extend(["", "## Automated observations", ""])
        report.extend(f"- {item}" for item in generated_observations)
    if plots:
        report.extend(["", "## Plots", ""])
        report.extend(f"- [{path.name}]({path.name})" for path in plots)
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    return report_path, summary_path, plots
