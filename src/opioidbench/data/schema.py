"""Typed schemas for raw, validated, and normalized EI-MS records."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import date
from typing import Any, Literal


@dataclass(frozen=True, slots=True)
class SpectrumPeak:
    """One centroid EI-MS peak."""

    mz: float
    intensity: float


@dataclass(frozen=True, slots=True)
class SourceProvenance:
    """Traceability fields attached to every spectrum."""

    source: str
    source_version: str | None
    source_record_id: str | None
    retrieval_date: date | None
    source_url: str | None
    license: str | None
    source_file: str | None = None
    source_file_sha256: str | None = None


@dataclass(frozen=True, slots=True)
class RawSpectrum:
    """Parsed spectrum that preserves source peak values and metadata."""

    spectrum_id: str
    compound_name_raw: str | None
    peaks: tuple[SpectrumPeak, ...]
    provenance: SourceProvenance
    metadata: dict[str, str] = field(default_factory=dict)
    analytical_modality: str = "ei_ms"
    ionization_method: str = "EI"


@dataclass(frozen=True, slots=True)
class ValidationIssue:
    """A machine-readable QC issue."""

    code: str
    severity: Literal["warning", "error"]
    message: str
    spectrum_id: str | None = None
    details: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class ValidatedSpectrum:
    """Raw spectrum paired with QC status; raw data remain untouched."""

    raw: RawSpectrum
    quality_flags: tuple[str, ...]
    is_valid: bool


@dataclass(frozen=True, slots=True)
class NormalizedSpectrumRecord:
    """Normalized, serializable research record for downstream curation."""

    spectrum_id: str
    compound_id: str
    source: str
    source_version: str | None
    source_record_id: str | None
    retrieval_date: str | None
    source_url: str | None
    license: str | None
    compound_name_raw: str | None
    compound_name_companion: str | None
    compound_name_resolution: str | None
    compound_name_canonical: str | None
    new_in_source_version: str | None
    cas: str | None
    pubchem_cid: str | None
    smiles: str | None
    canonical_smiles: str | None
    inchikey: str | None
    structure_source: str | None
    structure_source_record_id: str | None
    structure_match_method: str | None
    structure_formula_relation: str | None
    molecular_formula: str | None
    molecular_weight: float | None
    source_monoisotopic_mass: float | None
    source_m_plus_proton_mass: float | None
    analytical_modality: str
    ionization_method: str
    instrument_metadata: dict[str, str]
    opioid_label: str | None
    opioid_family: str | None
    scaffold_id: str | None
    murcko_scaffold_smiles: str | None
    quality_flags: tuple[str, ...]
    peaks: tuple[SpectrumPeak, ...]

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable dictionary."""
        return asdict(self)
