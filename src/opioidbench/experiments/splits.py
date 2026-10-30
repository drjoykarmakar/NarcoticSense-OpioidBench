"""Deterministic leakage-aware benchmark split generation and auditing."""

from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np


class BenchmarkDependencyError(RuntimeError):
    """Raised when optional benchmark dependencies are unavailable."""


@dataclass(frozen=True, slots=True)
class BenchmarkSample:
    spectrum_id: str
    compound_id: str
    label: int
    label_name: str
    scaffold_id: str | None
    structure_conflict: bool
    peaks: tuple[tuple[float, float], ...]
    analogue_group: str | None = None


def _operational_opioid_analogue_group(
    preferred_name: str,
    opioid_family: str,
    subfamily: str,
) -> str:
    """Map curated opioids to broad, auditable analogue lineages for holdout CV.

    This is an evaluation grouping, not an independent pharmacological label. The
    mapping deliberately uses broad lineages so related opioid chemistry cannot be
    split across train/test merely because exact Murcko scaffolds differ.
    """
    name = preferred_name.casefold().replace("α", "alpha")  # noqa: RUF001
    family = opioid_family.casefold()
    sub = subfamily.casefold()

    fentanyl_core_names = {"fentanyl", "alfentanil", "sufentanil", "remifentanil"}
    if "fentanyl" in name or name.strip() in fentanyl_core_names:
        return "fentanyl_anilidopiperidine"
    if "nitazene" in name or "nitazene" in sub:
        return "nitazene_benzimidazole"
    if family in {
        "natural opium alkaloids",
        "morphinan derivatives",
        "oripavine derivatives",
    }:
        return "morphinan_opiate"
    if family == "diphenylpropylamine derivatives":
        return "diphenylpropylamine"
    if family == "benzomorphan derivatives":
        return "benzomorphan"
    if family == "phenylpiperidine derivatives" or name.strip() == "anileridine":
        return "pethidine_phenylpiperidine"
    if name.strip() in {"tramadol", "tapentadol"}:
        return "tramadol_tapentadol"
    if name.strip() == "mt-45":
        return "mt45_piperazine"
    if "ap-237" in name:
        return "ap237_cinnamylpiperazine"
    if name.strip() == "u-47700":
        return "u_series_cyclohexylbenzamide"
    if name.strip() == "brorphine":
        return "brorphine_benzimidazolone"
    if name.strip() == "oliceridine":
        return "oliceridine"
    if name.strip() == "meptazinol":
        return "meptazinol"
    # Curated positives that do not fit one of the explicit broader lineages remain
    # isolated by name rather than being merged by an unsupported heuristic.
    normalized = "_".join(preferred_name.casefold().split())
    return f"singleton:{normalized}"


class _DisjointSet:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))
        self.rank = [0] * size

    def find(self, value: int) -> int:
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, left: int, right: int) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return
        if self.rank[left_root] < self.rank[right_root]:
            left_root, right_root = right_root, left_root
        self.parent[right_root] = left_root
        if self.rank[left_root] == self.rank[right_root]:
            self.rank[left_root] += 1


def _analogue_component_groups(samples: list[BenchmarkSample]) -> np.ndarray:
    """Union exact scaffolds and curated opioid analogue groups into holdout units."""
    dsu = _DisjointSet(len(samples))
    first_by_scaffold: dict[str, int] = {}
    first_by_analogue: dict[str, int] = {}
    for idx, sample in enumerate(samples):
        if sample.scaffold_id:
            first = first_by_scaffold.setdefault(sample.scaffold_id, idx)
            dsu.union(idx, first)
        if sample.label == 1 and sample.analogue_group:
            first = first_by_analogue.setdefault(sample.analogue_group, idx)
            dsu.union(idx, first)
    return np.asarray([f"component:{dsu.find(idx)}" for idx in range(len(samples))], dtype=object)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def load_curated_samples(
    dataset_path: Path,
    compounds_path: Path,
    taxonomy_path: Path,
    *,
    require_scaffold: bool = False,
) -> tuple[list[BenchmarkSample], dict[str, int | str]]:
    """Join normalized spectra to curated taxonomy and resolved compound metadata."""
    compounds = {
        row["compound_id"]: row for row in _load_csv(compounds_path) if row.get("compound_id")
    }
    taxonomy = {
        row["compound_id"]: row for row in _load_csv(taxonomy_path) if row.get("compound_id")
    }

    samples: list[BenchmarkSample] = []
    excluded_uncurated = 0
    excluded_missing_scaffold = 0
    excluded_structure_conflict = 0
    with dataset_path.open(encoding="utf-8") as handle:
        for line in handle:
            record = json.loads(line)
            compound_id = record.get("compound_id") or ""
            tax = taxonomy.get(compound_id)
            if not tax or tax.get("curation_status") != "curated":
                excluded_uncurated += 1
                continue
            label_name = tax.get("opioid_status") or ""
            if label_name not in {"opioid", "non-opioid"}:
                excluded_uncurated += 1
                continue
            compound = compounds.get(compound_id, {})
            conflict = (compound.get("structure_conflict") or "").strip().lower() == "true"
            scaffold = (compound.get("scaffold_id") or "").strip() or None
            if require_scaffold and conflict:
                excluded_structure_conflict += 1
                continue
            if require_scaffold and not scaffold:
                excluded_missing_scaffold += 1
                continue
            peaks = tuple(
                (float(peak["mz"]), float(peak["intensity"])) for peak in record.get("peaks", [])
            )
            analogue_group = None
            if label_name == "opioid":
                analogue_group = _operational_opioid_analogue_group(
                    tax.get("preferred_name") or "",
                    tax.get("opioid_family") or "",
                    tax.get("subfamily") or "",
                )
            samples.append(
                BenchmarkSample(
                    spectrum_id=str(record.get("spectrum_id") or ""),
                    compound_id=compound_id,
                    label=1 if label_name == "opioid" else 0,
                    label_name=label_name,
                    scaffold_id=scaffold,
                    structure_conflict=conflict,
                    peaks=peaks,
                    analogue_group=analogue_group,
                )
            )

    audit: dict[str, int | str] = {
        "dataset_sha256": _sha256_file(dataset_path),
        "compounds_sha256": _sha256_file(compounds_path),
        "taxonomy_sha256": _sha256_file(taxonomy_path),
        "included_spectra": len(samples),
        "included_opioid_spectra": sum(sample.label for sample in samples),
        "included_nonopioid_spectra": sum(1 - sample.label for sample in samples),
        "included_compounds": len({sample.compound_id for sample in samples}),
        "included_scaffolds": len({sample.scaffold_id for sample in samples if sample.scaffold_id}),
        "included_opioid_analogue_groups": len(
            {sample.analogue_group for sample in samples if sample.analogue_group}
        ),
        "excluded_uncurated_spectra": excluded_uncurated,
        "excluded_missing_scaffold_spectra": excluded_missing_scaffold,
        "excluded_structure_conflict_spectra": excluded_structure_conflict,
    }
    return samples, audit


def _spectral_matrix(samples: list[BenchmarkSample], mz_max: int = 700) -> np.ndarray:
    matrix = np.zeros((len(samples), mz_max + 1), dtype=np.float32)
    for row_idx, sample in enumerate(samples):
        for mz, intensity in sample.peaks:
            bin_idx = round(mz)
            if 0 <= bin_idx <= mz_max:
                matrix[row_idx, bin_idx] += float(intensity)
        norm = float(np.linalg.norm(matrix[row_idx]))
        if norm:
            matrix[row_idx] /= norm
    return matrix


def _peak_hash(sample: BenchmarkSample) -> str:
    payload = json.dumps(sample.peaks, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _fold_leakage_audit(
    samples: list[BenchmarkSample],
    fold_assignment: np.ndarray,
    *,
    near_identical_threshold: float = 0.9999,
) -> list[dict[str, int | float]]:
    matrix = _spectral_matrix(samples)
    peak_hashes = [_peak_hash(sample) for sample in samples]
    reports: list[dict[str, int | float]] = []
    for fold in sorted(set(int(value) for value in fold_assignment.tolist())):
        test_idx = np.flatnonzero(fold_assignment == fold)
        train_idx = np.flatnonzero(fold_assignment != fold)
        train_compounds = {samples[idx].compound_id for idx in train_idx}
        test_compounds = {samples[idx].compound_id for idx in test_idx}
        train_scaffolds = {
            samples[idx].scaffold_id for idx in train_idx if samples[idx].scaffold_id
        }
        test_scaffolds = {samples[idx].scaffold_id for idx in test_idx if samples[idx].scaffold_id}
        train_analogue_groups = {
            samples[idx].analogue_group
            for idx in train_idx
            if samples[idx].label == 1 and samples[idx].analogue_group
        }
        test_analogue_groups = {
            samples[idx].analogue_group
            for idx in test_idx
            if samples[idx].label == 1 and samples[idx].analogue_group
        }
        train_hashes = {peak_hashes[idx] for idx in train_idx}
        test_hashes = {peak_hashes[idx] for idx in test_idx}

        near_pairs = 0
        if len(test_idx) and len(train_idx):
            similarities = matrix[test_idx] @ matrix[train_idx].T
            near_pairs = int(np.count_nonzero(similarities >= near_identical_threshold))

        reports.append(
            {
                "fold": fold,
                "train_spectra": len(train_idx),
                "test_spectra": len(test_idx),
                "train_opioids": int(sum(samples[idx].label for idx in train_idx)),
                "test_opioids": int(sum(samples[idx].label for idx in test_idx)),
                "train_nonopioids": int(sum(1 - samples[idx].label for idx in train_idx)),
                "test_nonopioids": int(sum(1 - samples[idx].label for idx in test_idx)),
                "compound_overlap": len(train_compounds & test_compounds),
                "scaffold_overlap": len(train_scaffolds & test_scaffolds),
                "opioid_analogue_group_overlap": len(train_analogue_groups & test_analogue_groups),
                "exact_spectrum_overlap": len(train_hashes & test_hashes),
                "near_identical_cross_split_pairs": near_pairs,
                "near_identical_threshold": near_identical_threshold,
            }
        )
    return reports


def generate_cv_splits(
    dataset_path: Path,
    compounds_path: Path,
    taxonomy_path: Path,
    *,
    output_dir: Path,
    protocol: str,
    n_splits: int = 5,
    seed: int = 20260929,
    near_identical_threshold: float = 0.9999,
) -> dict[str, object]:
    """Generate frozen CV fold assignments and aggregate leakage audits."""
    try:
        from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold
    except ImportError as exc:  # pragma: no cover - dependency dependent
        raise BenchmarkDependencyError(
            "scikit-learn is required for benchmark split generation; install opioidbench[models]."
        ) from exc

    if protocol not in {"random", "compound", "scaffold", "analogue"}:
        raise ValueError("protocol must be one of: random, compound, scaffold, analogue")
    if n_splits < 2:
        raise ValueError("n_splits must be >= 2")

    samples, input_audit = load_curated_samples(
        dataset_path,
        compounds_path,
        taxonomy_path,
        require_scaffold=protocol in {"scaffold", "analogue"},
    )
    if not samples:
        raise ValueError("No curated samples are available for split generation")
    y = np.asarray([sample.label for sample in samples], dtype=np.int8)
    indices = np.arange(len(samples))
    assignments = np.full(len(samples), -1, dtype=np.int16)

    if protocol == "random":
        splitter = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        iterator = splitter.split(indices, y)
    else:
        if protocol == "compound":
            groups = np.asarray([sample.compound_id for sample in samples], dtype=object)
        elif protocol == "scaffold":
            groups = np.asarray([sample.scaffold_id for sample in samples], dtype=object)
        else:
            groups = _analogue_component_groups(samples)
        splitter = StratifiedGroupKFold(
            n_splits=n_splits,
            shuffle=True,
            random_state=seed,
        )
        iterator = splitter.split(indices, y, groups)

    for fold, (_, test_idx) in enumerate(iterator):
        assignments[test_idx] = fold
    if np.any(assignments < 0):
        raise RuntimeError("Split assignment left one or more samples unassigned")

    leakage = _fold_leakage_audit(
        samples,
        assignments,
        near_identical_threshold=near_identical_threshold,
    )
    if protocol in {"compound", "scaffold", "analogue"} and any(
        int(item["compound_overlap"]) for item in leakage
    ):
        raise RuntimeError("Compound leakage detected in a group-held-out protocol")
    if protocol in {"scaffold", "analogue"} and any(
        int(item["scaffold_overlap"]) for item in leakage
    ):
        raise RuntimeError("Scaffold leakage detected in scaffold/analogue-held-out protocol")
    if protocol == "analogue" and any(
        int(item["opioid_analogue_group_overlap"]) for item in leakage
    ):
        raise RuntimeError("Opioid analogue-group leakage detected in analogue-held-out protocol")

    output_dir.mkdir(parents=True, exist_ok=True)
    assignments_path = output_dir / f"{protocol}_cv{n_splits}_assignments.csv"
    with assignments_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "spectrum_id",
                "compound_id",
                "label",
                "scaffold_id",
                "analogue_group",
                "fold",
            ],
            lineterminator="\n",
        )
        writer.writeheader()
        for sample, fold in zip(samples, assignments, strict=True):
            writer.writerow(
                {
                    "spectrum_id": sample.spectrum_id,
                    "compound_id": sample.compound_id,
                    "label": sample.label_name,
                    "scaffold_id": sample.scaffold_id or "",
                    "analogue_group": sample.analogue_group or "",
                    "fold": int(fold),
                }
            )

    report: dict[str, object] = {
        "protocol": protocol,
        "n_splits": n_splits,
        "seed": seed,
        "near_identical_threshold": near_identical_threshold,
        "input": input_audit,
        "folds": leakage,
        "assignments_file": assignments_path.name,
        "assignments_sha256": _sha256_file(assignments_path),
        "interpretation": (
            "Frozen pilot CV on curated labels only. Unknown/pending compounds are excluded. "
            "Random CV is diagnostic; compound/scaffold/analogue CV are leakage-aware. "
            "Analogue CV unions exact scaffolds with broad curated opioid analogue lineages."
        ),
    }
    report_path = output_dir / f"{protocol}_cv{n_splits}_audit.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def load_fold_assignments(path: Path) -> dict[str, int]:
    """Load spectrum_id -> fold assignments produced by :func:`generate_cv_splits`."""
    rows = _load_csv(path)
    return {row["spectrum_id"]: int(row["fold"]) for row in rows}
