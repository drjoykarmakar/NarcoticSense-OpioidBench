"""Exploratory nonlinear EI-MS models and train-reference similarity diagnostics.

This is an additive, unfrozen research experiment. No threshold, model selection,
or tuning is performed using the held-out test fold.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from opioidbench.experiments.baselines import METRIC_NAMES, _metrics, vectorize_eims
from opioidbench.experiments.splits import BenchmarkDependencyError, BenchmarkSample


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def evaluate_advanced(
    samples: list[BenchmarkSample],
    fold_by_spectrum: dict[str, int],
    *,
    output_dir: Path,
    mz_max: int = 700,
    seed: int = 20260929,
    input_hashes: dict[str, str] | None = None,
) -> dict[str, object]:
    """Run predeclared RBF SVM and ExtraTrees on frozen held-out folds.

    Per-sample predictions are written only to the caller's local output directory.
    The maximum train-spectrum cosine is an uncalibrated novelty *proxy*, not
    validated out-of-distribution detection or a calibrated probability.
    """
    try:
        from sklearn.ensemble import ExtraTreesClassifier
        from sklearn.svm import SVC
    except ImportError as exc:  # pragma: no cover
        raise BenchmarkDependencyError(
            "Install opioidbench[models] for advanced model evaluation."
        ) from exc

    if not samples:
        raise ValueError("No curated samples supplied")
    ids = [sample.spectrum_id for sample in samples]
    if len(set(ids)) != len(ids):
        raise ValueError("Duplicate spectrum IDs are not allowed")
    if any(sample.label not in {0, 1} for sample in samples):
        raise ValueError("Only curated binary labels are supported")
    if any(sid not in fold_by_spectrum for sid in ids):
        raise ValueError("Every sample requires a frozen fold assignment")
    if mz_max < 1:
        raise ValueError("mz_max must be positive")

    folds = np.asarray([fold_by_spectrum[sid] for sid in ids], dtype=int)
    fold_ids = sorted(set(folds.tolist()))
    if len(fold_ids) < 2:
        raise ValueError("At least two folds are required")
    x = vectorize_eims(samples, mz_max=mz_max, sqrt_transform=True)
    y = np.asarray([sample.label for sample in samples], dtype=np.int8)

    fold_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    for fold in fold_ids:
        train = np.flatnonzero(folds != fold)
        test = np.flatnonzero(folds == fold)
        if len(np.unique(y[train])) != 2 or len(np.unique(y[test])) != 2:
            raise ValueError(f"Fold {fold} must contain both classes in train and test")

        x_train, x_test = x[train], x[test]
        y_train, y_test = y[train], y[test]
        # Max cosine to *training* references only. No held-out data enters fit.
        nearest_train_cosine = np.max(x_test @ x_train.T, axis=1)

        models = {
            "rbf_svm": SVC(
                C=1.0, kernel="rbf", gamma="scale",
                class_weight="balanced", random_state=seed,
            ),
            "extra_trees": ExtraTreesClassifier(
                n_estimators=500, max_features="sqrt", min_samples_leaf=1,
                class_weight="balanced", random_state=seed, n_jobs=1,
            ),
        }
        for name, model in models.items():
            model.fit(x_train, y_train)
            y_pred = model.predict(x_test).astype(np.int8)
            if name == "rbf_svm":
                scores = np.asarray(model.decision_function(x_test), dtype=float)
            else:
                scores = np.asarray(model.predict_proba(x_test)[:, 1], dtype=float)
            metrics = _metrics(y_test, y_pred, scores)
            fold_rows.append({
                "model": name, "fold": fold, "n_train": len(train),
                "n_test": len(test), **metrics,
                "mean_nearest_train_cosine": float(np.mean(nearest_train_cosine)),
            })
            for i, sample_idx in enumerate(test):
                prediction_rows.append({
                    "model": name, "fold": fold,
                    "spectrum_id": samples[sample_idx].spectrum_id,
                    "compound_id": samples[sample_idx].compound_id,
                    "true_label": samples[sample_idx].label_name,
                    "predicted_label": "opioid" if y_pred[i] else "non-opioid",
                    "score": float(scores[i]),
                    "nearest_train_cosine": float(nearest_train_cosine[i]),
                })

    aggregate: dict[str, dict[str, float]] = {}
    for name in ("rbf_svm", "extra_trees"):
        rows = [row for row in fold_rows if row["model"] == name]
        summary: dict[str, float] = {}
        for metric in (*METRIC_NAMES, "mean_nearest_train_cosine"):
            values = np.asarray([row[metric] for row in rows], dtype=float)
            summary[f"{metric}_mean"] = float(np.nanmean(values))
            summary[f"{metric}_std"] = float(np.nanstd(values, ddof=0))
        aggregate[name] = summary

    output_dir.mkdir(parents=True, exist_ok=True)
    fold_fields = [
        "model", "fold", "n_train", "n_test", *METRIC_NAMES,
        "tn", "fp", "fn", "tp", "mean_nearest_train_cosine",
    ]
    with (output_dir / "fold_metrics.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fold_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(fold_rows)

    prediction_fields = [
        "model", "fold", "spectrum_id", "compound_id", "true_label",
        "predicted_label", "score", "nearest_train_cosine",
    ]
    with (output_dir / "predictions.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=prediction_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(prediction_rows)

    report: dict[str, object] = {
        "status": "exploratory_unfrozen",
        "n_samples": len(samples),
        "n_opioid": int(np.sum(y)),
        "n_nonopioid": int(len(y) - np.sum(y)),
        "n_folds": len(fold_ids),
        "mz_max": mz_max,
        "seed": seed,
        "input_sha256": input_hashes or {},
        "preprocessing": "nominal integer m/z; sqrt intensity; L2 normalization",
        "models": {
            "rbf_svm": "SVC RBF; C=1; gamma=scale; class_weight=balanced",
            "extra_trees": "500 trees; sqrt features; class_weight=balanced",
        },
        "novelty_diagnostic": (
            "Nearest cosine similarity to the training fold only. This is not "
            "calibrated uncertainty, validated OOD detection, or an abstention rule."
        ),
        "aggregate": aggregate,
        "fold_metrics_file": "fold_metrics.csv",
        "predictions_file": "predictions.csv",
        "publication_warning": (
            "Exploratory results; not a replacement for the frozen v0.4.0 "
            "benchmark. Resolve scaffold-provenance discrepancy before comparison."
        ),
    }
    (output_dir / "summary.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return report
