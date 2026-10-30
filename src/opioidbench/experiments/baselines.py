"""Simple, fixed-configuration EI-MS baselines for pilot benchmarking."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np

from opioidbench.experiments.splits import BenchmarkDependencyError, BenchmarkSample

METRIC_NAMES = (
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


def vectorize_eims(
    samples: list[BenchmarkSample],
    *,
    mz_max: int = 700,
    sqrt_transform: bool = True,
) -> np.ndarray:
    """Bin nominal m/z, optionally square-root intensities, and L2-normalize."""
    matrix = np.zeros((len(samples), mz_max + 1), dtype=np.float32)
    for row_idx, sample in enumerate(samples):
        for mz, intensity in sample.peaks:
            bin_idx = round(mz)
            if 0 <= bin_idx <= mz_max:
                matrix[row_idx, bin_idx] += max(0.0, float(intensity))
        if sqrt_transform:
            np.sqrt(matrix[row_idx], out=matrix[row_idx])
        norm = float(np.linalg.norm(matrix[row_idx]))
        if norm:
            matrix[row_idx] /= norm
    return matrix


def _metrics(y_true: np.ndarray, y_pred: np.ndarray, score: np.ndarray) -> dict[str, float]:
    try:
        from sklearn.metrics import (
            accuracy_score,
            average_precision_score,
            balanced_accuracy_score,
            f1_score,
            matthews_corrcoef,
            precision_score,
            recall_score,
            roc_auc_score,
        )
    except ImportError as exc:  # pragma: no cover
        raise BenchmarkDependencyError(
            "scikit-learn is required for benchmark evaluation; install opioidbench[models]."
        ) from exc

    tn = int(np.sum((y_true == 0) & (y_pred == 0)))
    fp = int(np.sum((y_true == 0) & (y_pred == 1)))
    specificity = tn / (tn + fp) if (tn + fp) else math.nan
    result = {
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "sensitivity": float(recall_score(y_true, y_pred, zero_division=0)),
        "specificity": float(specificity),
        "precision": float(precision_score(y_true, y_pred, zero_division=0)),
        "f1": float(f1_score(y_true, y_pred, zero_division=0)),
        "mcc": float(matthews_corrcoef(y_true, y_pred)),
        "auroc": (
            float(roc_auc_score(y_true, score)) if len(set(y_true.tolist())) == 2 else math.nan
        ),
        "average_precision": (
            float(average_precision_score(y_true, score))
            if len(set(y_true.tolist())) == 2
            else math.nan
        ),
        "tn": tn,
        "fp": fp,
        "fn": int(np.sum((y_true == 1) & (y_pred == 0))),
        "tp": int(np.sum((y_true == 1) & (y_pred == 1))),
    }
    return result


def _nearest_neighbor_scores(
    x_train: np.ndarray,
    y_train: np.ndarray,
    x_test: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    similarities = x_test @ x_train.T
    positive = np.flatnonzero(y_train == 1)
    negative = np.flatnonzero(y_train == 0)
    if not len(positive) or not len(negative):
        raise ValueError("Nearest-neighbor baseline requires both classes in the training fold")
    pos_score = np.max(similarities[:, positive], axis=1)
    neg_score = np.max(similarities[:, negative], axis=1)
    decision = pos_score - neg_score
    return (decision > 0).astype(np.int8), decision


def evaluate_baselines(
    samples: list[BenchmarkSample],
    fold_by_spectrum: dict[str, int],
    *,
    output_dir: Path,
    mz_max: int = 700,
    seed: int = 20260929,
) -> dict[str, object]:
    """Evaluate fixed baselines across a frozen fold assignment."""
    try:
        from sklearn.ensemble import RandomForestClassifier
        from sklearn.linear_model import LogisticRegression
        from sklearn.svm import LinearSVC
    except ImportError as exc:  # pragma: no cover
        raise BenchmarkDependencyError(
            "scikit-learn is required for benchmark evaluation; install opioidbench[models]."
        ) from exc

    missing = [
        sample.spectrum_id for sample in samples if sample.spectrum_id not in fold_by_spectrum
    ]
    if missing:
        raise ValueError(f"Fold assignment is missing {len(missing)} curated spectrum(s)")

    x = vectorize_eims(samples, mz_max=mz_max, sqrt_transform=True)
    y = np.asarray([sample.label for sample in samples], dtype=np.int8)
    folds = np.asarray([fold_by_spectrum[sample.spectrum_id] for sample in samples], dtype=np.int16)
    unique_folds = sorted(set(int(value) for value in folds.tolist()))

    fold_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    model_names = ("spectral_1nn", "logistic_regression", "linear_svm", "random_forest")

    for fold in unique_folds:
        train_idx = np.flatnonzero(folds != fold)
        test_idx = np.flatnonzero(folds == fold)
        x_train, y_train = x[train_idx], y[train_idx]
        x_test, y_test = x[test_idx], y[test_idx]
        if len(set(y_train.tolist())) != 2 or len(set(y_test.tolist())) != 2:
            raise ValueError(f"Fold {fold} does not contain both classes in train and test")

        predictions: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        predictions["spectral_1nn"] = _nearest_neighbor_scores(x_train, y_train, x_test)

        logreg = LogisticRegression(
            C=1.0,
            class_weight="balanced",
            max_iter=5000,
            solver="liblinear",
            random_state=seed,
        )
        logreg.fit(x_train, y_train)
        predictions["logistic_regression"] = (
            logreg.predict(x_test).astype(np.int8),
            logreg.decision_function(x_test).astype(float),
        )

        svm = LinearSVC(C=1.0, class_weight="balanced", random_state=seed)
        svm.fit(x_train, y_train)
        predictions["linear_svm"] = (
            svm.predict(x_test).astype(np.int8),
            svm.decision_function(x_test).astype(float),
        )

        forest = RandomForestClassifier(
            n_estimators=500,
            max_features="sqrt",
            min_samples_leaf=1,
            class_weight="balanced_subsample",
            random_state=seed,
            n_jobs=1,
        )
        forest.fit(x_train, y_train)
        predictions["random_forest"] = (
            forest.predict(x_test).astype(np.int8),
            forest.predict_proba(x_test)[:, 1].astype(float),
        )

        for model_name in model_names:
            y_pred, score = predictions[model_name]
            metrics = _metrics(y_test, y_pred, score)
            fold_rows.append(
                {
                    "model": model_name,
                    "fold": fold,
                    "n_train": len(train_idx),
                    "n_test": len(test_idx),
                    **metrics,
                }
            )
            for local_idx, sample_idx in enumerate(test_idx):
                prediction_rows.append(
                    {
                        "model": model_name,
                        "fold": fold,
                        "spectrum_id": samples[sample_idx].spectrum_id,
                        "compound_id": samples[sample_idx].compound_id,
                        "true_label": samples[sample_idx].label_name,
                        "predicted_label": "opioid" if int(y_pred[local_idx]) else "non-opioid",
                        "score": float(score[local_idx]),
                    }
                )

    aggregate: dict[str, dict[str, float]] = {}
    for model_name in model_names:
        model_rows = [row for row in fold_rows if row["model"] == model_name]
        summary: dict[str, float] = {}
        for metric in METRIC_NAMES:
            values = np.asarray([float(row[metric]) for row in model_rows], dtype=float)
            summary[f"{metric}_mean"] = float(np.nanmean(values))
            summary[f"{metric}_std"] = float(np.nanstd(values, ddof=0))
        aggregate[model_name] = summary

    output_dir.mkdir(parents=True, exist_ok=True)
    fold_path = output_dir / "fold_metrics.csv"
    fieldnames = ["model", "fold", "n_train", "n_test", *METRIC_NAMES, "tn", "fp", "fn", "tp"]
    with fold_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(fold_rows)

    prediction_path = output_dir / "predictions.csv"
    with prediction_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "model",
                "fold",
                "spectrum_id",
                "compound_id",
                "true_label",
                "predicted_label",
                "score",
            ],
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(prediction_rows)

    report: dict[str, object] = {
        "n_samples": len(samples),
        "n_opioid": int(np.sum(y)),
        "n_nonopioid": int(len(y) - np.sum(y)),
        "n_folds": len(unique_folds),
        "mz_max": mz_max,
        "preprocessing": "nominal integer m/z bins; sqrt intensity; L2 normalization",
        "models": {
            "spectral_1nn": "max cosine to each class; decision=max_opioid-max_nonopioid",
            "logistic_regression": "C=1, class_weight=balanced, liblinear",
            "linear_svm": "LinearSVC C=1, class_weight=balanced",
            "random_forest": "500 trees, sqrt features, class_weight=balanced_subsample",
        },
        "aggregate": aggregate,
        "fold_metrics_file": fold_path.name,
        "predictions_file": prediction_path.name,
        "interpretation": (
            "Pilot baseline on the evidence-curated subset only. These metrics are not the final "
            "benchmark and must not be generalized to the 3,801-compound source inventory."
        ),
    }
    (output_dir / "summary.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return report
