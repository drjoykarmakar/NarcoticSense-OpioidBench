"""Cross-source EI-MS overlap auditing without redistributing source spectra."""

from __future__ import annotations

import hashlib
import json
import math
import zipfile
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from opioidbench.data.enfsi import (
    _conservative_aliases,
    _load_normalized_records,
    _normalize_name,
    read_enfsi_entry_list,
)
from opioidbench.data.parsers import parse_jcamp_file, parse_jcamp_text
from opioidbench.data.schema import RawSpectrum
from opioidbench.data.validation import validate_dataset

DEFAULT_OVERLAP_THRESHOLDS = (0.99999, 0.9999, 0.999, 0.99, 0.95)


@dataclass(frozen=True, slots=True)
class SpectralMatchDetail:
    """Local-only cross-source match detail containing no peak list or compound name."""

    enfsi_source_record_id: str | None
    swgdrug_source_record_id: str | None
    cosine_similarity: float
    mz_support_jaccard: float
    identical_mz_support: bool
    exact_peak_list: bool
    near_identical: bool
    metadata_overlap_status: str | None
    metadata_best_match_consistent: bool | None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class ENFSISpectralOverlapReport:
    """Aggregate-only ENFSI-versus-SWGDRUG spectrum-overlap report."""

    enfsi_spectra: int
    enfsi_valid_spectra: int
    enfsi_invalid_spectra: int
    enfsi_duplicate_spectrum_groups: int
    enfsi_duplicate_spectrum_records: int
    enfsi_placeholder_cas_records: int
    enfsi_invalid_cas_records: int
    enfsi_invalid_molecular_weight_records: int
    swgdrug_spectra: int
    exact_peak_list_match_records: int
    identical_mz_support_match_records: int
    near_identical_cosine_threshold: float
    near_identical_best_match_records: int
    spectrum_disjoint_candidate_records: int
    best_match_threshold_counts: dict[str, int]
    best_match_cosine_quantiles: dict[str, float]
    entry_list_alignment_verified: bool | None
    entry_list_record_id_matches: int | None
    entry_list_formula_matches: int | None
    entry_list_molecular_weight_matches: int | None
    high_confidence_metadata_matched_entries: int | None
    ambiguous_metadata_entries: int | None
    unmatched_metadata_entries: int | None
    near_identical_high_confidence_metadata_matched: int | None
    near_identical_ambiguous_metadata: int | None
    near_identical_metadata_unmatched: int | None
    metadata_unmatched_spectrum_disjoint_candidates: int | None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _load_numpy():
    try:
        import numpy as np
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise RuntimeError(
            "Spectrum-overlap auditing requires NumPy. Install opioidbench[spectral]."
        ) from exc
    return np


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def read_jcamp_source(path: Path, *, source: str, source_version: str) -> tuple[
    list[RawSpectrum], dict[str, str | None]
]:
    """Read a JCAMP/JDX file or a single-member ZIP without persisting restricted data."""
    if not path.is_file():
        raise FileNotFoundError(path)
    if path.suffix.lower() != ".zip":
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        return (
            parse_jcamp_file(
                path,
                source=source,
                source_version=source_version,
                source_file_sha256=digest,
            ),
            {
                "archive_sha256": None,
                "member": path.name,
                "member_sha256": digest,
            },
        )

    archive_payload = path.read_bytes()
    archive_digest = _sha256_bytes(archive_payload)
    allowed = {".jdx", ".dx", ".hpj"}
    with zipfile.ZipFile(path) as archive:
        candidates = [
            info
            for info in archive.infolist()
            if not info.is_dir() and Path(info.filename).suffix.lower() in allowed
        ]
        if len(candidates) != 1:
            names = [info.filename for info in candidates]
            raise ValueError(
                "Spectrum archive must contain exactly one JDX/DX/HPJ member; "
                f"found {len(candidates)}: {names}"
            )
        info = candidates[0]
        member_payload = archive.read(info)
    member_digest = _sha256_bytes(member_payload)
    try:
        text = member_payload.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ValueError(f"Spectrum member {info.filename!r} is not UTF-8 text") from exc
    records = parse_jcamp_text(
        text,
        source=source,
        source_version=source_version,
        source_file=info.filename,
        source_file_sha256=member_digest,
    )
    return records, {
        "archive_sha256": archive_digest,
        "member": info.filename,
        "member_sha256": member_digest,
    }


def _peak_signature(record: RawSpectrum) -> tuple[tuple[float, float], ...]:
    return tuple((peak.mz, peak.intensity) for peak in record.peaks)


def _mz_signature(record: RawSpectrum) -> tuple[float, ...]:
    return tuple(peak.mz for peak in record.peaks)


def _support_jaccard(first: RawSpectrum, second: RawSpectrum) -> float:
    first_mz = {peak.mz for peak in first.peaks}
    second_mz = {peak.mz for peak in second.peaks}
    union = first_mz | second_mz
    if not union:
        return 1.0
    return len(first_mz & second_mz) / len(union)


def _cosine_best_matches(
    queries: list[RawSpectrum],
    references: list[RawSpectrum],
    *,
    chunk_size: int = 256,
) -> tuple[list[int], list[float]]:
    """Return the nearest reference and raw-intensity cosine for each query.

    Intensities are L2-normalized, so source-specific global intensity scales such as
    999 versus 9,999 do not affect the score. Nominal integer m/z bins are used because
    both audited EI-MS libraries contain centroided integer m/z values.
    """
    if not queries or not references:
        raise ValueError("Spectrum-overlap audit requires non-empty query and reference sets")
    np = _load_numpy()
    max_mz = int(
        max(
            max(peak.mz for peak in record.peaks)
            for record in [*queries, *references]
            if record.peaks
        )
    )

    def vectorize(records: list[RawSpectrum]):
        matrix = np.zeros((len(records), max_mz + 1), dtype=np.float64)
        for row, record in enumerate(records):
            for peak in record.peaks:
                mz = int(peak.mz)
                if float(mz) != peak.mz:
                    raise ValueError(
                        "Cross-source overlap audit currently requires integer nominal m/z values"
                    )
                matrix[row, mz] = float(peak.intensity)
        norms = np.linalg.norm(matrix, axis=1, keepdims=True)
        if bool(np.any(norms == 0)):
            raise ValueError("Cannot compare an empty or all-zero spectrum")
        return matrix / norms

    reference_matrix = vectorize(references)
    best_indices: list[int] = []
    best_scores: list[float] = []
    for start in range(0, len(queries), chunk_size):
        chunk = vectorize(queries[start : start + chunk_size])
        similarities = chunk @ reference_matrix.T
        indices = np.argmax(similarities, axis=1)
        scores = similarities[np.arange(similarities.shape[0]), indices]
        best_indices.extend(int(value) for value in indices.tolist())
        best_scores.extend(float(value) for value in scores.tolist())
    return best_indices, best_scores


def _quantiles(values: list[float]) -> dict[str, float]:
    np = _load_numpy()
    points = (0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.0)
    return {
        f"q{int(point * 100):02d}": round(float(np.quantile(values, point)), 9)
        for point in points
    }


def _metadata_match_candidates(
    entry_list_path: Path,
    normalized_records_path: Path,
) -> tuple[list[set[str]], list[dict[str, Any]]]:
    entries = read_enfsi_entry_list(entry_list_path)
    normalized = _load_normalized_records(normalized_records_path)
    index: dict[tuple[str, str | None, int | None], set[str]] = defaultdict(set)
    for record in normalized:
        weight = record.get("molecular_weight")
        rounded_weight = round(float(weight)) if weight is not None else None
        names = {
            record.get("compound_name_raw"),
            record.get("compound_name_companion"),
            record.get("compound_name_canonical"),
        }
        for name in names:
            if isinstance(name, str) and name.strip():
                index[
                    (_normalize_name(name), record.get("molecular_formula"), rounded_weight)
                ].add(str(record.get("compound_id") or record.get("spectrum_id") or ""))

    candidates: list[set[str]] = []
    for entry in entries:
        matched: set[str] = set()
        rounded_weight = round(entry.molecular_weight)
        for alias in _conservative_aliases(entry.name):
            matched.update(
                index.get((_normalize_name(alias), entry.formula, rounded_weight), set())
            )
        candidates.append(matched)
    return candidates, normalized


def _verify_entry_alignment(
    spectra: list[RawSpectrum],
    entry_list_path: Path,
) -> tuple[bool, int, int, int]:
    entries = read_enfsi_entry_list(entry_list_path)
    if len(entries) != len(spectra):
        return False, 0, 0, 0
    record_id_matches = 0
    formula_matches = 0
    molecular_weight_matches = 0
    for spectrum, entry in zip(spectra, entries, strict=True):
        if spectrum.provenance.source_record_id == f"Library Entry {entry.source_index}":
            record_id_matches += 1
        if spectrum.metadata.get("MOLFORM") == entry.formula:
            formula_matches += 1
        try:
            weight = float(spectrum.metadata.get("MW", ""))
        except ValueError:
            weight = math.nan
        if weight == entry.molecular_weight:
            molecular_weight_matches += 1
    verified = record_id_matches == len(entries) and formula_matches == len(entries)
    return verified, record_id_matches, formula_matches, molecular_weight_matches


def audit_enfsi_spectral_overlap(
    enfsi_spectra: list[RawSpectrum],
    swgdrug_spectra: list[RawSpectrum],
    *,
    near_identical_threshold: float = 0.9999,
    entry_list_path: Path | None = None,
    normalized_swgdrug_path: Path | None = None,
) -> tuple[ENFSISpectralOverlapReport, list[SpectralMatchDetail]]:
    """Audit cross-source spectral overlap and return aggregate and local-only details.

    A ``near_identical`` flag is a leakage-risk screen, not a claim that two measurements
    are literally the same acquisition. The default 0.9999 threshold is intentionally
    stringent and is accompanied by a threshold-sensitivity table in the report.
    """
    if not 0 < near_identical_threshold <= 1:
        raise ValueError("near-identical threshold must lie in (0, 1]")
    _, enfsi_qc = validate_dataset(enfsi_spectra)
    if enfsi_qc.number_invalid:
        raise ValueError(
            f"ENFSI overlap audit requires valid spectra; found {enfsi_qc.number_invalid} invalid"
        )

    exact_reference_signatures = {_peak_signature(record) for record in swgdrug_spectra}
    mz_reference_signatures = {_mz_signature(record) for record in swgdrug_spectra}
    exact_flags = [
        bool(record.peaks) and _peak_signature(record) in exact_reference_signatures
        for record in enfsi_spectra
    ]
    mz_flags = [
        bool(record.peaks) and _mz_signature(record) in mz_reference_signatures
        for record in enfsi_spectra
    ]

    best_indices, best_scores = _cosine_best_matches(enfsi_spectra, swgdrug_spectra)
    near_flags = [score >= near_identical_threshold for score in best_scores]
    threshold_counts = {
        f"{threshold:.5f}".rstrip("0").rstrip("."): sum(
            score >= threshold for score in best_scores
        )
        for threshold in DEFAULT_OVERLAP_THRESHOLDS
    }

    entry_verified: bool | None = None
    entry_record_id_matches: int | None = None
    entry_formula_matches: int | None = None
    entry_weight_matches: int | None = None
    metadata_candidates: list[set[str]] | None = None
    normalized_records: list[dict[str, Any]] | None = None
    if entry_list_path is not None:
        (
            entry_verified,
            entry_record_id_matches,
            entry_formula_matches,
            entry_weight_matches,
        ) = _verify_entry_alignment(enfsi_spectra, entry_list_path)
        if not entry_verified:
            raise ValueError(
                "ENFSI JDX does not align safely with the supplied entry list by record ID/formula"
            )
    if (entry_list_path is None) != (normalized_swgdrug_path is None):
        raise ValueError(
            "entry_list_path and normalized_swgdrug_path must be supplied together for "
            "joint metadata/spectrum auditing"
        )
    if entry_list_path is not None and normalized_swgdrug_path is not None:
        metadata_candidates, normalized_records = _metadata_match_candidates(
            entry_list_path, normalized_swgdrug_path
        )
        if len(metadata_candidates) != len(enfsi_spectra):
            raise ValueError("ENFSI entry list and JDX record counts differ")
        if len(normalized_records) != len(swgdrug_spectra):
            raise ValueError("Normalized SWGDRUG records and raw SWGDRUG spectra differ in count")

    details: list[SpectralMatchDetail] = []
    metadata_status_counts: Counter[str] = Counter()
    near_metadata_status_counts: Counter[str] = Counter()
    unmatched_spectrum_disjoint = 0
    for index, (best_index, score) in enumerate(zip(best_indices, best_scores, strict=True)):
        status: str | None = None
        consistent: bool | None = None
        if metadata_candidates is not None and normalized_records is not None:
            candidates = metadata_candidates[index]
            if len(candidates) == 1:
                status = "unique"
            elif len(candidates) > 1:
                status = "ambiguous"
            else:
                status = "unmatched"
            metadata_status_counts[status] += 1
            if near_flags[index]:
                near_metadata_status_counts[status] += 1
            if status == "unmatched" and not near_flags[index]:
                unmatched_spectrum_disjoint += 1
            if status == "unique":
                consistent = str(normalized_records[best_index].get("compound_id")) in candidates

        details.append(
            SpectralMatchDetail(
                enfsi_source_record_id=enfsi_spectra[index].provenance.source_record_id,
                swgdrug_source_record_id=swgdrug_spectra[best_index].provenance.source_record_id,
                cosine_similarity=round(score, 12),
                mz_support_jaccard=round(
                    _support_jaccard(enfsi_spectra[index], swgdrug_spectra[best_index]), 12
                ),
                identical_mz_support=mz_flags[index],
                exact_peak_list=exact_flags[index],
                near_identical=near_flags[index],
                metadata_overlap_status=status,
                metadata_best_match_consistent=consistent,
            )
        )

    duplicate_spectrum_records = sum(
        len(issue.details.get("spectrum_ids", []))
        for issue in enfsi_qc.issues
        if issue.code == "duplicate_spectrum"
    )
    report = ENFSISpectralOverlapReport(
        enfsi_spectra=len(enfsi_spectra),
        enfsi_valid_spectra=enfsi_qc.number_valid,
        enfsi_invalid_spectra=enfsi_qc.number_invalid,
        enfsi_duplicate_spectrum_groups=enfsi_qc.issue_counts.get("duplicate_spectrum", 0),
        enfsi_duplicate_spectrum_records=duplicate_spectrum_records,
        enfsi_placeholder_cas_records=enfsi_qc.issue_counts.get("placeholder_cas", 0),
        enfsi_invalid_cas_records=enfsi_qc.issue_counts.get("invalid_cas", 0),
        enfsi_invalid_molecular_weight_records=enfsi_qc.issue_counts.get(
            "invalid_molecular_weight", 0
        ),
        swgdrug_spectra=len(swgdrug_spectra),
        exact_peak_list_match_records=sum(exact_flags),
        identical_mz_support_match_records=sum(mz_flags),
        near_identical_cosine_threshold=near_identical_threshold,
        near_identical_best_match_records=sum(near_flags),
        spectrum_disjoint_candidate_records=sum(not flag for flag in near_flags),
        best_match_threshold_counts=threshold_counts,
        best_match_cosine_quantiles=_quantiles(best_scores),
        entry_list_alignment_verified=entry_verified,
        entry_list_record_id_matches=entry_record_id_matches,
        entry_list_formula_matches=entry_formula_matches,
        entry_list_molecular_weight_matches=entry_weight_matches,
        high_confidence_metadata_matched_entries=(
            metadata_status_counts.get("unique", 0) if metadata_candidates is not None else None
        ),
        ambiguous_metadata_entries=(
            metadata_status_counts.get("ambiguous", 0)
            if metadata_candidates is not None
            else None
        ),
        unmatched_metadata_entries=(
            metadata_status_counts.get("unmatched", 0)
            if metadata_candidates is not None
            else None
        ),
        near_identical_high_confidence_metadata_matched=(
            near_metadata_status_counts.get("unique", 0)
            if metadata_candidates is not None
            else None
        ),
        near_identical_ambiguous_metadata=(
            near_metadata_status_counts.get("ambiguous", 0)
            if metadata_candidates is not None
            else None
        ),
        near_identical_metadata_unmatched=(
            near_metadata_status_counts.get("unmatched", 0)
            if metadata_candidates is not None
            else None
        ),
        metadata_unmatched_spectrum_disjoint_candidates=(
            unmatched_spectrum_disjoint if metadata_candidates is not None else None
        ),
    )
    return report, details


def write_private_match_manifest(path: Path, details: list[SpectralMatchDetail]) -> None:
    """Write a local-only match manifest containing IDs/scores but no spectra or names."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for detail in details:
            handle.write(json.dumps(detail.to_dict(), sort_keys=True) + "\n")
