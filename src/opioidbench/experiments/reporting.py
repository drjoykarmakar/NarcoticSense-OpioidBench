"""Validate local pilot outputs and export aggregate-only benchmark reports.

Never reads predictions.csv, assignments, spectra, or compound identifiers.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

METRICS = (
    "accuracy",
    "balanced_accuracy",
    "sensitivity",
    "specificity",
    "precision",
    "f1",
    "mcc",
    "auroc",
    "average_precision",
)
MODELS = ("spectral_1nn", "logistic_regression", "linear_svm", "random_forest")
PROTOCOLS = ("random", "compound", "scaffold", "analogue")
COUNTS = ("tn", "fp", "fn", "tp")
TOLERANCE = 1e-8


def _number(value: object, label: str) -> float:
    try:
        result = float(value)
    except (ValueError, TypeError) as exc:
        raise ValueError(f"Invalid numeric value for {label}") from exc
    if not math.isfinite(result):
        raise ValueError(f"Non-finite value for {label}")
    return result


def _integer(value: object, label: str) -> int:
    number = _number(value, label)
    if not number.is_integer() or number < 0:
        raise ValueError(f"Invalid nonnegative integer for {label}")
    return int(number)


def _close(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=TOLERANCE, abs_tol=TOLERANCE)


def _read_run(protocol: str, directory: Path) -> dict[str, object]:
    if protocol not in PROTOCOLS:
        raise ValueError(f"Unsupported protocol: {protocol}")
    summary_path = directory / "summary.json"
    metrics_path = directory / "fold_metrics.csv"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if not isinstance(summary, dict) or not isinstance(summary.get("aggregate"), dict):
        raise ValueError(f"Malformed summary.json for {protocol}")

    with metrics_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        expected = {"model", "fold", "n_train", "n_test", *METRICS, *COUNTS}
        if reader.fieldnames is None or not expected.issubset(reader.fieldnames):
            raise ValueError(f"Missing required fold metric columns for {protocol}")
        rows = list(reader)

    if not rows:
        raise ValueError(f"No fold metrics for {protocol}")
    n_folds = _integer(summary.get("n_folds"), "n_folds")
    n_samples = _integer(summary.get("n_samples"), "n_samples")
    n_opioid = _integer(summary.get("n_opioid"), "n_opioid")
    n_nonopioid = _integer(summary.get("n_nonopioid"), "n_nonopioid")
    if n_folds < 2 or n_samples == 0 or n_opioid + n_nonopioid != n_samples:
        raise ValueError(f"Inconsistent cohort sizes for {protocol}")
    if set(summary["aggregate"]) != set(MODELS):
        raise ValueError(f"Missing or unexpected models for {protocol}")

    output_models: dict[str, object] = {}
    for model in MODELS:
        model_rows = [row for row in rows if row["model"] == model]
        if len(model_rows) != n_folds:
            raise ValueError(f"Expected {n_folds} folds for {protocol}/{model}")
        folds = [_integer(row["fold"], "fold") for row in model_rows]
        if len(set(folds)) != n_folds:
            raise ValueError(f"Duplicate folds for {protocol}/{model}")
        if set(folds) != set(range(n_folds)):
            raise ValueError(f"Unexpected fold indices for {protocol}/{model}")
        n_test_total = 0
        totals = {key: 0 for key in COUNTS}
        values: dict[str, list[float]] = {metric: [] for metric in METRICS}
        for row in model_rows:
            n_test = _integer(row["n_test"], "n_test")
            n_train = _integer(row["n_train"], "n_train")
            if n_test == 0 or n_train + n_test != n_samples:
                raise ValueError(f"Invalid fold size for {protocol}/{model}")
            n_test_total += n_test
            counts = {key: _integer(row[key], key) for key in COUNTS}
            if sum(counts.values()) != n_test:
                raise ValueError(f"Confusion counts do not match test size: {protocol}/{model}")
            for key in COUNTS:
                totals[key] += counts[key]
            for metric in METRICS:
                val = _number(row[metric], metric)
                if metric != "mcc" and not 0 <= val <= 1:
                    raise ValueError(f"Metric out of range: {protocol}/{model}/{metric}")
                if metric == "mcc" and not -1 <= val <= 1:
                    raise ValueError(f"MCC out of range: {protocol}/{model}")
                values[metric].append(val)
        if n_test_total != n_samples:
            raise ValueError(f"Test folds do not partition cohort for {protocol}/{model}")
        if totals["tp"] + totals["fn"] != n_opioid:
            raise ValueError(f"Opioid totals inconsistent for {protocol}/{model}")
        if totals["tn"] + totals["fp"] != n_nonopioid:
            raise ValueError(f"Comparator totals inconsistent for {protocol}/{model}")
        means = {}
        stds = {}
        for metric, series in values.items():
            mean = sum(series) / len(series)
            std = math.sqrt(sum((value - mean) ** 2 for value in series) / len(series))
            aggregate = summary["aggregate"][model]
            if not _close(mean, _number(aggregate.get(f"{metric}_mean"), metric)):
                raise ValueError(f"Mean mismatch: {protocol}/{model}/{metric}")
            if not _close(std, _number(aggregate.get(f"{metric}_std"), metric)):
                raise ValueError(f"Std mismatch: {protocol}/{model}/{metric}")
            means[metric] = mean
            stds[metric] = std
        output_models[model] = {
            "metrics_mean": means,
            "metrics_std": stds,
            "confusion_counts_pooled": totals,
        }
    if len(rows) != len(MODELS) * n_folds:
        raise ValueError(f"Unexpected fold rows for {protocol}")
    return {
        "n_samples": n_samples,
        "n_opioid": n_opioid,
        "n_nonopioid": n_nonopioid,
        "n_folds": n_folds,
        "models": output_models,
    }


def generate_public_report(runs: dict[str, Path], output_dir: Path) -> dict[str, object]:
    """Generate a minimal allowlisted aggregate report from local experiment directories.

    Inputs may contain private files, but only summary.json and fold_metrics.csv
    are opened. Output paths should be reviewed before publication.
    """
    if not runs:
        raise ValueError("At least one protocol run is required")
    if any(protocol not in PROTOCOLS for protocol in runs):
        raise ValueError("Unknown protocol")
    reports = {protocol: _read_run(protocol, Path(path)) for protocol, path in sorted(runs.items())}
    public = {
        "schema_version": "1.0",
        "status": "preliminary evidence-curated pilot; not externally validated",
        "protocols": reports,
        "metric_aggregation": "unweighted mean and population std over held-out folds",
        "pooled_confusion_note": (
            "Confusion counts sum out-of-fold test predictions; "
            "pooled metrics may differ from fold-mean metrics."
        ),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "pilot_metrics.json").write_text(
        json.dumps(public, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    fields = ["protocol", "model", "n_samples", "n_opioid", "n_nonopioid", "n_folds"]
    fields += [f"{metric}_{stat}" for metric in METRICS for stat in ("mean", "std")]
    fields += list(COUNTS)
    with (output_dir / "pilot_metrics.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for protocol, report in reports.items():
            for model, model_report in report["models"].items():
                row = {"protocol": protocol, "model": model}
                row.update({key: report[key] for key in fields if key in report})
                for metric in METRICS:
                    row[f"{metric}_mean"] = model_report["metrics_mean"][metric]
                    row[f"{metric}_std"] = model_report["metrics_std"][metric]
                row.update(model_report["confusion_counts_pooled"])
                writer.writerow(row)
    return public
