"""Synthetic-only tests for the public aggregate reporting pipeline."""

import csv
import json

import pytest

from opioidbench.experiments.reporting import METRICS, MODELS, generate_public_report


def make_run(path, *, corrupt=False):
    path.mkdir()
    fields = ["model", "fold", "n_train", "n_test", *METRICS, "tn", "fp", "fn", "tp"]
    rows = []
    aggregate = {}
    for model in MODELS:
        for fold in range(2):
            row = {
                "model": model,
                "fold": fold,
                "n_train": 4,
                "n_test": 4,
                "tn": 2,
                "fp": 0,
                "fn": 0,
                "tp": 2,
            }
            row.update({metric: 1.0 for metric in METRICS})
            rows.append(row)
        aggregate[model] = {
            f"{metric}_{stat}": (1.0 if stat == "mean" else 0.0)
            for metric in METRICS
            for stat in ("mean", "std")
        }
    if corrupt:
        rows[0]["tp"] = 3
    with (path / "fold_metrics.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    (path / "summary.json").write_text(
        json.dumps(
            {
                "n_folds": 2,
                "n_samples": 8,
                "n_opioid": 4,
                "n_nonopioid": 4,
                "aggregate": aggregate,
                "predictions_file": "secret_predictions.csv",
                "private_identifier": "must_not_escape",
            }
        ),
        encoding="utf-8",
    )
    (path / "predictions.csv").write_text("SECRET_COMPOUND_ID,restricted_spectrum\n")
    return path


def test_report_is_aggregate_only(tmp_path):
    run = make_run(tmp_path / "local")
    result = generate_public_report({"compound": run}, tmp_path / "results" / "report")
    assert result["protocols"]["compound"]["models"]["random_forest"][
        "confusion_counts_pooled"
    ] == {"tn": 4, "fp": 0, "fn": 0, "tp": 4}
    output = tmp_path / "results" / "report"
    assert (output / "pilot_metrics.csv").is_file()
    exported = (output / "pilot_metrics.json").read_text()
    assert "SECRET_COMPOUND_ID" not in exported
    assert "must_not_escape" not in exported
    assert "secret_predictions.csv" not in exported
    assert "compound_id" not in exported


def test_report_rejects_inconsistent_confusion_counts(tmp_path):
    run = make_run(tmp_path / "local", corrupt=True)
    with pytest.raises(ValueError, match="Confusion counts"):
        generate_public_report({"compound": run}, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_report_rejects_unknown_protocol(tmp_path):
    with pytest.raises(ValueError, match="Unknown protocol"):
        generate_public_report({"external": tmp_path}, tmp_path / "out")


def test_report_rejects_missing_fold(tmp_path):
    run = make_run(tmp_path / "local")
    path = run / "fold_metrics.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    fields = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows[:-1])
    with pytest.raises(ValueError, match="Expected 2 folds"):
        generate_public_report({"compound": run}, tmp_path / "out")


def test_report_rejects_tampered_summary(tmp_path):
    run = make_run(tmp_path / "local")
    path = run / "summary.json"
    data = json.loads(path.read_text())
    data["aggregate"]["random_forest"]["f1_mean"] = 0.3
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="Mean mismatch"):
        generate_public_report({"compound": run}, tmp_path / "out")
