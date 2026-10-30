"""Automated QC for raw EI-MS records."""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from typing import Any, Literal

from opioidbench.chemistry.identifiers import (
    CAS_PLACEHOLDERS,
    cas_checksum_valid,
    conservative_compound_id,
    molecular_formula_plausible,
    normalize_cas,
    normalize_compound_name,
)
from opioidbench.data.companion import preferred_source_name
from opioidbench.data.schema import (
    NormalizedSpectrumRecord,
    RawSpectrum,
    ValidatedSpectrum,
    ValidationIssue,
)


@dataclass(frozen=True, slots=True)
class QCReport:
    """Machine-readable dataset QC report."""

    number_spectra: int
    number_valid: int
    number_invalid: int
    issue_counts: dict[str, int]
    issues: tuple[ValidationIssue, ...]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def validate_spectrum(spectrum: RawSpectrum) -> tuple[ValidatedSpectrum, list[ValidationIssue]]:
    """Validate one raw spectrum without mutating it."""
    issues: list[ValidationIssue] = []

    def add(
        code: str,
        severity: Literal["warning", "error"],
        message: str,
        **details: Any,
    ) -> None:
        issues.append(
            ValidationIssue(
                code=code,
                severity=severity,
                message=message,
                spectrum_id=spectrum.spectrum_id,
                details=details,
            )
        )

    if not spectrum.compound_name_raw or not spectrum.compound_name_raw.strip():
        add("missing_compound_name", "error", "Compound name is missing")
    if not spectrum.peaks:
        add("empty_spectrum", "error", "Spectrum contains no peaks")

    declared_points = spectrum.metadata.get("NPOINTS")
    if declared_points:
        try:
            expected_points = int(declared_points.strip())
        except ValueError:
            add(
                "invalid_npoints",
                "warning",
                "NPOINTS is not an integer",
                value=declared_points,
            )
        else:
            if expected_points != len(spectrum.peaks):
                add(
                    "npoints_mismatch",
                    "error",
                    "Declared NPOINTS does not match parsed peak count",
                    declared=expected_points,
                    parsed=len(spectrum.peaks),
                )

    raw_cas = spectrum.metadata.get("CAS_REGISTRY_NO") or spectrum.metadata.get("CAS")
    if raw_cas:
        if raw_cas.strip() in CAS_PLACEHOLDERS:
            add(
                "placeholder_cas",
                "warning",
                "Source uses a placeholder CAS value; normalized CAS remains missing",
                value=raw_cas,
            )
        elif not cas_checksum_valid(raw_cas):
            add(
                "invalid_cas",
                "warning",
                "Source CAS value fails format/check-digit validation",
                value=raw_cas,
            )

    formula = spectrum.metadata.get("MOLFORM") or spectrum.metadata.get("FORMULA")
    if formula and not molecular_formula_plausible(formula):
        add(
            "suspicious_molecular_formula",
            "warning",
            "Molecular formula has unsupported or suspicious syntax",
            value=formula,
        )

    molecular_weight = None
    source_mw = spectrum.metadata.get("MW")
    if source_mw:
        try:
            molecular_weight = float(source_mw)
        except ValueError:
            add(
                "invalid_molecular_weight",
                "warning",
                "Source molecular weight is not numeric",
                value=source_mw,
            )
        else:
            if molecular_weight <= 0:
                add(
                    "invalid_molecular_weight",
                    "warning",
                    "Source molecular weight must be positive",
                    value=source_mw,
                )

    source_exact_mass = spectrum.metadata.get("EXACT_MASS_LIST_MONOISOTOPIC_MASS")
    if source_exact_mass:
        try:
            exact_mass = float(source_exact_mass)
        except ValueError:
            add(
                "invalid_source_monoisotopic_mass",
                "warning",
                "Source monoisotopic mass is not numeric",
                value=source_exact_mass,
            )
        else:
            if exact_mass <= 0:
                add(
                    "invalid_source_monoisotopic_mass",
                    "warning",
                    "Source monoisotopic mass must be positive",
                    value=source_exact_mass,
                )
            elif molecular_weight is not None and round(exact_mass) != molecular_weight:
                add(
                    "exact_mass_nominal_weight_rounding_mismatch",
                    "warning",
                    "Rounded source monoisotopic mass differs from source nominal molecular weight",
                    monoisotopic_mass=exact_mass,
                    molecular_weight=molecular_weight,
                )

    mz_counts: Counter[float] = Counter()
    for index, peak in enumerate(spectrum.peaks):
        if peak.mz <= 0:
            add("invalid_mz", "error", "m/z must be positive", index=index, mz=peak.mz)
        if peak.intensity < 0:
            add(
                "negative_intensity",
                "error",
                "Intensity must be non-negative",
                index=index,
                intensity=peak.intensity,
            )
        mz_counts[peak.mz] += 1

    duplicate_mz = sorted(mz for mz, count in mz_counts.items() if count > 1)
    if duplicate_mz:
        add(
            "duplicate_peaks",
            "warning",
            "Spectrum contains repeated m/z values",
            mz_values=duplicate_mz,
        )

    flags = tuple(sorted({issue.code for issue in issues}))
    is_valid = not any(issue.severity == "error" for issue in issues)
    return ValidatedSpectrum(raw=spectrum, quality_flags=flags, is_valid=is_valid), issues


def validate_dataset(spectra: list[RawSpectrum]) -> tuple[list[ValidatedSpectrum], QCReport]:
    """Validate records and detect dataset-level duplicates/conflicts."""
    validated: list[ValidatedSpectrum] = []
    all_issues: list[ValidationIssue] = []
    for spectrum in spectra:
        item, issues = validate_spectrum(spectrum)
        validated.append(item)
        all_issues.extend(issues)

    id_groups: dict[str, list[RawSpectrum]] = defaultdict(list)
    signature_groups: dict[tuple[tuple[float, float], ...], list[RawSpectrum]] = defaultdict(list)
    source_record_groups: dict[tuple[str, str | None], list[RawSpectrum]] = defaultdict(list)
    compound_groups: dict[str, list[RawSpectrum]] = defaultdict(list)
    for spectrum in spectra:
        id_groups[spectrum.spectrum_id].append(spectrum)
        signature = tuple((p.mz, p.intensity) for p in spectrum.peaks)
        signature_groups[signature].append(spectrum)
        source_key = (spectrum.provenance.source, spectrum.provenance.source_record_id)
        source_record_groups[source_key].append(spectrum)
        if spectrum.compound_name_raw:
            compound_id = conservative_compound_id(
                cas=spectrum.metadata.get("CAS_REGISTRY_NO") or spectrum.metadata.get("CAS"),
                name=preferred_source_name(spectrum),
            )
            compound_groups[compound_id].append(spectrum)

    dataset_issues: list[ValidationIssue] = []
    for spectrum_id, group in id_groups.items():
        if len(group) > 1:
            dataset_issues.append(
                ValidationIssue(
                    code="duplicate_spectrum_id",
                    severity="error",
                    message="Duplicate spectrum_id occurs in the dataset",
                    spectrum_id=spectrum_id,
                    details={"count": len(group)},
                )
            )
    for signature, group in signature_groups.items():
        if signature and len(group) > 1:
            dataset_issues.append(
                ValidationIssue(
                    code="duplicate_spectrum",
                    severity="warning",
                    message="Identical peak lists occur in multiple records",
                    details={"spectrum_ids": [item.spectrum_id for item in group]},
                )
            )
    for (source, record_id), group in source_record_groups.items():
        if record_id and len(group) > 1:
            dataset_issues.append(
                ValidationIssue(
                    code="duplicate_source_record_id",
                    severity="warning",
                    message="A source record identifier occurs more than once",
                    details={
                        "source": source,
                        "source_record_id": record_id,
                        "spectrum_ids": [item.spectrum_id for item in group],
                    },
                )
            )
    for compound_id, group in compound_groups.items():
        if len(group) > 1:
            dataset_issues.append(
                ValidationIssue(
                    code="duplicate_chemical_identity",
                    severity="warning",
                    message="Multiple spectra resolve to the same conservative chemical identity",
                    details={
                        "compound_id": compound_id,
                        "count": len(group),
                        "spectrum_ids": [item.spectrum_id for item in group],
                    },
                )
            )
        names = {
            normalize_compound_name(preferred_source_name(item))
            for item in group
            if normalize_compound_name(preferred_source_name(item))
        }
        if len(names) > 1:
            dataset_issues.append(
                ValidationIssue(
                    code="conflicting_compound_names",
                    severity="warning",
                    message="One conservative chemical identity has multiple source names",
                    details={
                        "compound_id": compound_id,
                        "names": sorted(names),
                        "source_record_ids": [item.provenance.source_record_id for item in group],
                    },
                )
            )

    all_issues.extend(dataset_issues)
    invalid_ids = {
        issue.spectrum_id for issue in all_issues if issue.severity == "error" and issue.spectrum_id
    }
    # Rebuild flags so dataset-level errors are reflected in record status.
    rebuilt: list[ValidatedSpectrum] = []
    per_spectrum_codes: dict[str, set[str]] = defaultdict(set)
    for issue in all_issues:
        if issue.spectrum_id:
            per_spectrum_codes[issue.spectrum_id].add(issue.code)
    for item in validated:
        sid = item.raw.spectrum_id
        rebuilt.append(
            ValidatedSpectrum(
                raw=item.raw,
                quality_flags=tuple(sorted(per_spectrum_codes[sid])),
                is_valid=item.is_valid and sid not in invalid_ids,
            )
        )

    counts = Counter(issue.code for issue in all_issues)
    report = QCReport(
        number_spectra=len(spectra),
        number_valid=sum(item.is_valid for item in rebuilt),
        number_invalid=sum(not item.is_valid for item in rebuilt),
        issue_counts=dict(sorted(counts.items())),
        issues=tuple(all_issues),
    )
    return rebuilt, report


def _optional_positive_float(value: str | None) -> float | None:
    if value is None:
        return None
    try:
        parsed = float(value)
    except ValueError:
        return None
    return parsed if parsed > 0 else None


def normalize_spectrum(spectrum: ValidatedSpectrum) -> NormalizedSpectrumRecord:
    """Create the milestone normalized record without fabricating chemical metadata."""
    raw = spectrum.raw
    metadata = raw.metadata
    raw_cas = metadata.get("CAS_REGISTRY_NO") or metadata.get("CAS")
    companion_name = metadata.get("COMPOUND_LIST_NAME")
    name_resolution = metadata.get("COMPOUND_LIST_NAME_RELATION")
    source_name = preferred_source_name(raw)
    canonical_name = normalize_compound_name(source_name)
    compound_id = conservative_compound_id(cas=raw_cas, name=canonical_name)
    source_mw = metadata.get("MW")
    molecular_weight = None
    if source_mw:
        try:
            parsed_mw = float(source_mw)
        except ValueError:
            pass
        else:
            if parsed_mw > 0:
                molecular_weight = parsed_mw
    return NormalizedSpectrumRecord(
        spectrum_id=raw.spectrum_id,
        compound_id=compound_id,
        source=raw.provenance.source,
        source_version=raw.provenance.source_version,
        source_record_id=raw.provenance.source_record_id,
        retrieval_date=(
            raw.provenance.retrieval_date.isoformat() if raw.provenance.retrieval_date else None
        ),
        source_url=raw.provenance.source_url,
        license=raw.provenance.license,
        compound_name_raw=raw.compound_name_raw,
        compound_name_companion=companion_name,
        compound_name_resolution=name_resolution,
        compound_name_canonical=canonical_name,
        new_in_source_version=metadata.get("COMPOUND_LIST_NEW_IN_VERSION"),
        cas=normalize_cas(raw_cas),
        pubchem_cid=None,
        smiles=metadata.get("NIST_STRUCTURE_CANONICAL_SMILES") or None,
        canonical_smiles=metadata.get("NIST_STRUCTURE_CANONICAL_SMILES") or None,
        inchikey=metadata.get("NIST_STRUCTURE_INCHIKEY") or None,
        structure_source=(
            "swgdrug_nist_usrstruc" if metadata.get("NIST_STRUCTURE_CANONICAL_SMILES") else None
        ),
        structure_source_record_id=metadata.get("NIST_STRUCTURE_LIBRARY_ID") or None,
        structure_match_method=metadata.get("NIST_STRUCTURE_MATCH_METHOD") or None,
        structure_formula_relation=metadata.get("NIST_STRUCTURE_FORMULA_RELATION") or None,
        molecular_formula=metadata.get("MOLFORM") or metadata.get("FORMULA"),
        molecular_weight=molecular_weight,
        source_monoisotopic_mass=_optional_positive_float(
            metadata.get("EXACT_MASS_LIST_MONOISOTOPIC_MASS")
        ),
        source_m_plus_proton_mass=_optional_positive_float(
            metadata.get("EXACT_MASS_LIST_M_PLUS_PROTON")
        ),
        analytical_modality=raw.analytical_modality,
        ionization_method=raw.ionization_method,
        instrument_metadata={
            key: value
            for key, value in metadata.items()
            if key
            in {
                "ORIGIN",
                "OWNER",
                "SPECTROMETER_DATA_SYSTEM",
                "INSTRUMENT",
                "SAMPLE_DESCRIPTION",
                "SOURCE_COMMENT_FULL",
                "SOURCE_COMMENT_RELATION",
                "SOURCE_COMMENT_SOURCE_FORMAT",
            }
        },
        opioid_label=None,
        opioid_family=None,
        scaffold_id=metadata.get("NIST_STRUCTURE_SCAFFOLD_ID") or None,
        murcko_scaffold_smiles=(
            metadata.get("NIST_STRUCTURE_SCAFFOLD_SMILES")
            if metadata.get("NIST_STRUCTURE_CANONICAL_SMILES")
            else None
        ),
        quality_flags=spectrum.quality_flags,
        peaks=raw.peaks,
    )
