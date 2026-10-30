"""Tests for exploratory nonlinear models (synthetic data only)."""

import csv
import json
from pathlib import Path

import pytest

from opioidbench.experiments.advanced import evaluate_advanced
from opioidbench.experiments.splits import BenchmarkSample


def _samples():
    samples = []
    folds = {}
    for fold in range(5):
        for label in (0, 1):
            for repeat in range(2):
                sid = f"synthetic_{fold}_{label}_{repeat}"
                samples.append(
                    BenchmarkSample(
                        spectrum_id=sid,
                        compound_id=sid,
                        label=label,
                        label_name="opioid" if label else "non-opioid",
                        scaffold_id=f"scaffold_{fold}_{label}_{repeat}",
                        structure_conflict=False,
                        peaks=((50.0 if label == 0 else 200.0, 100.0),),
                    )
                )
                folds[sid] = fold
    return samples, folds


def test_advanced_models_run_and_report_train_only_similarity(tmp_path: Path):
    samples, folds = _samples()
    report = evaluate_advanced(
        samples, folds, output_dir=tmp_path / "advanced", mz_max=700, seed=20260929
    )
    assert report["status"] == "exploratory_unfrozen"
    assert report["n_samples"] == 20
    assert set(report["aggregate"]) == {"rbf_svm", "extra_trees"}
    assert (tmp_path / "advanced" / "summary.json").is_file()
    with (tmp_path / "advanced" / "predictions.csv").open(newline="", encoding="utf-8") as f:
        predictions = list(csv.DictReader(f))
    assert len(predictions) == 40
    assert all(0.0 <= float(row["nearest_train_cosine"]) <= 1.000001 for row in predictions)
    assert all(row["spectrum_id"].startswith("synthetic_") for row in predictions)
    with (tmp_path / "advanced" / "fold_metrics.csv").open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 10
    assert json.loads((tmp_path / "advanced" / "summary.json").read_text())["mz_max"] == 700


def test_advanced_rejects_missing_folds(tmp_path: Path):
    samples, folds = _samples()
    del folds[samples[0].spectrum_id]
    with pytest.raises(ValueError, match="frozen fold assignment"):
        evaluate_advanced(samples, folds, output_dir=tmp_path)


def test_advanced_rejects_single_class_fold(tmp_path: Path):
    samples, folds = _samples()
    for sample in samples:
        folds[sample.spectrum_id] = sample.label
    with pytest.raises(ValueError, match="both classes"):
        evaluate_advanced(samples, folds, output_dir=tmp_path)


def test_advanced_rejects_duplicate_spectrum_ids(tmp_path: Path):
    samples, folds = _samples()
    with pytest.raises(ValueError, match="Duplicate spectrum"):
        evaluate_advanced([*samples, samples[0]], folds, output_dir=tmp_path)
