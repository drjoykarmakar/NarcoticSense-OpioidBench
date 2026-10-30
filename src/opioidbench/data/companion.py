"""Companion metadata alignment for SWGDRUG release files."""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass, replace
from pathlib import Path

from opioidbench.data.schema import RawSpectrum


@dataclass(frozen=True, slots=True)
class CompoundListEntry:
    """One source-provided compound-list row."""

    row_number: int
    name: str
    new_in_source_version: str | None


@dataclass(frozen=True, slots=True)
class CompanionAlignmentReport:
    """Machine-readable result of aligning a companion list to spectra by source order."""

    number_entries: int
    relation_counts: dict[str, int]
    number_new_in_source_version: int
    source_version_marker: str | None

    def to_dict(self) -> dict[str, object]:
        return {
            "number_entries": self.number_entries,
            "relation_counts": dict(self.relation_counts),
            "number_new_in_source_version": self.number_new_in_source_version,
            "source_version_marker": self.source_version_marker,
        }


@dataclass(frozen=True, slots=True)
class ExactMassEntry:
    """One row from the SWGDRUG formula/exact-mass workbook."""

    row_number: int
    formula: str
    name: str
    monoisotopic_mass: float
    m_plus_proton_mass: float | None


@dataclass(frozen=True, slots=True)
class ExactMassAlignmentReport:
    """Alignment summary for the source-provided exact-mass companion workbook."""

    number_entries: int
    name_relation_counts: dict[str, int]
    number_formula_matches: int
    number_m_plus_proton_values: int
    number_m_plus_proton_missing: int
    source_version_marker: str | None

    def to_dict(self) -> dict[str, object]:
        return {
            "number_entries": self.number_entries,
            "name_relation_counts": dict(self.name_relation_counts),
            "number_formula_matches": self.number_formula_matches,
            "number_m_plus_proton_values": self.number_m_plus_proton_values,
            "number_m_plus_proton_missing": self.number_m_plus_proton_missing,
            "source_version_marker": self.source_version_marker,
        }


def _load_openpyxl():
    try:
        from openpyxl import load_workbook
    except ImportError as exc:  # pragma: no cover - depends on optional environment
        raise RuntimeError(
            "Reading SWGDRUG XLSX companion files requires the optional 'excel' dependency. "
            "Install opioidbench[excel]."
        ) from exc
    return load_workbook


def _version_marker(header: str) -> str | None:
    match = re.search(r"\bv\s*([0-9]+(?:\.[0-9]+)*)\b", header, flags=re.IGNORECASE)
    if match:
        return match.group(1)
    match = re.search(r"\bversion\s*([0-9]+(?:\.[0-9]+)*)\b", header, flags=re.IGNORECASE)
    return match.group(1) if match else None


def _is_red_font(cell: object) -> bool:
    color = getattr(getattr(cell, "font", None), "color", None)
    if color is None or getattr(color, "type", None) != "rgb":
        return False
    rgb = getattr(color, "rgb", None)
    return isinstance(rgb, str) and rgb.upper().endswith("FF0000")


def read_swgdrug_compound_list(path: Path) -> tuple[list[CompoundListEntry], str | None]:
    """Read the one-column SWGDRUG compound list while preserving red-entry semantics."""
    load_workbook = _load_openpyxl()
    workbook = load_workbook(path, read_only=False, data_only=True)
    if not workbook.sheetnames:
        raise ValueError(f"Workbook {path} contains no worksheets")
    sheet = workbook[workbook.sheetnames[0]]
    if sheet.max_column != 1:
        raise ValueError(
            f"Expected a one-column SWGDRUG compound list, found {sheet.max_column} columns"
        )
    header_value = sheet.cell(1, 1).value
    header = str(header_value).strip() if header_value is not None else ""
    if "compound name" not in header.casefold():
        raise ValueError(f"Unexpected SWGDRUG compound-list header: {header!r}")
    version = _version_marker(header)

    entries: list[CompoundListEntry] = []
    for row_number in range(2, sheet.max_row + 1):
        cell = sheet.cell(row_number, 1)
        if cell.value is None or not str(cell.value).strip():
            raise ValueError(f"Blank compound-list entry at row {row_number}")
        entries.append(
            CompoundListEntry(
                row_number=row_number,
                name=str(cell.value).strip(),
                new_in_source_version=version if _is_red_font(cell) else None,
            )
        )
    return entries, version


def read_swgdrug_exact_mass_list(path: Path) -> tuple[list[ExactMassEntry], str | None]:
    """Read SWGDRUG's source-distributed formula/exact-mass workbook.

    The workbook describes its masses as formula-derived by Agilent PCDL and unverified.
    This reader therefore preserves them as source metadata rather than treating them as
    independently validated chemical measurements.
    """
    load_workbook = _load_openpyxl()
    workbook = load_workbook(path, read_only=True, data_only=True)
    if not workbook.sheetnames:
        raise ValueError(f"Workbook {path} contains no worksheets")
    sheet = workbook[workbook.sheetnames[0]]

    header_row = None
    for row_number, row in enumerate(
        sheet.iter_rows(min_row=1, max_row=min(sheet.max_row, 20), max_col=4, values_only=True),
        start=1,
    ):
        values = [str(value).strip() if value is not None else "" for value in row]
        if values[:4] == ["Formula", "Compound Name", "Mass", "M+proton"]:
            header_row = row_number
            break
    if header_row is None:
        raise ValueError("Could not locate Formula/Compound Name/Mass/M+proton header row")

    version = None
    for row in sheet.iter_rows(min_row=1, max_row=header_row - 1, max_col=4, values_only=True):
        joined = " ".join(str(value) for value in row if value is not None)
        marker = _version_marker(joined)
        if marker:
            version = marker
            break

    entries: list[ExactMassEntry] = []
    for row_number, row in enumerate(
        sheet.iter_rows(min_row=header_row + 1, max_col=4, values_only=True),
        start=header_row + 1,
    ):
        formula, name, mass, m_plus_proton = row
        if all(value is None for value in row):
            continue
        if formula is None or not str(formula).strip():
            raise ValueError(f"Blank formula in exact-mass workbook row {row_number}")
        if name is None or not str(name).strip():
            raise ValueError(f"Blank compound name in exact-mass workbook row {row_number}")
        if mass is None:
            raise ValueError(f"Missing monoisotopic mass in exact-mass workbook row {row_number}")
        try:
            parsed_mass = float(mass)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Invalid monoisotopic mass in exact-mass workbook row {row_number}: {mass!r}"
            ) from exc
        if parsed_mass <= 0:
            raise ValueError(
                f"Non-positive monoisotopic mass in exact-mass workbook row {row_number}: {mass!r}"
            )
        parsed_proton = None
        if m_plus_proton is not None and str(m_plus_proton).strip():
            try:
                parsed_proton = float(m_plus_proton)
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"Invalid M+proton value in exact-mass workbook row {row_number}: "
                    f"{m_plus_proton!r}"
                ) from exc
        entries.append(
            ExactMassEntry(
                row_number=row_number,
                formula=str(formula).strip(),
                name=str(name).strip(),
                monoisotopic_mass=parsed_mass,
                m_plus_proton_mass=parsed_proton,
            )
        )
    return entries, version


def classify_name_relation(hpj_name: str, companion_name: str) -> str:
    """Classify how an HPJ name relates to its source companion-list name."""
    hpj = hpj_name.strip()
    companion = companion_name.strip()
    if hpj == companion:
        return "exact"
    if len(hpj) == 65 and companion.startswith(hpj):
        return "hpj_truncated"
    if hpj.startswith(companion):
        return "companion_shorter"
    return "mismatch"


def preferred_source_name(spectrum: RawSpectrum) -> str | None:
    """Return the source name used for identity fallback without altering raw source text."""
    companion = spectrum.metadata.get("COMPOUND_LIST_NAME")
    relation = spectrum.metadata.get("COMPOUND_LIST_NAME_RELATION")
    if companion and relation == "hpj_truncated":
        return companion
    return spectrum.compound_name_raw


def align_swgdrug_compound_list(
    spectra: list[RawSpectrum],
    entries: list[CompoundListEntry],
    *,
    strict: bool = True,
) -> tuple[list[RawSpectrum], CompanionAlignmentReport]:
    """Align companion rows to spectra in release order and attach auditable metadata.

    SWGDRUG v3.14's HPJ and companion workbook have the same release order. The
    function validates every row relationship before allowing the companion name to
    affect downstream identity fallback. Only the observed 65-character HPJ prefix
    pattern is treated as a high-confidence truncation repair.
    """
    if len(spectra) != len(entries):
        raise ValueError(
            "Spectrum/compound-list row count mismatch: "
            f"{len(spectra)} spectra vs {len(entries)} companion entries"
        )

    relations: list[str] = []
    for index, (spectrum, entry) in enumerate(zip(spectra, entries, strict=True), start=1):
        if not spectrum.compound_name_raw:
            relation = "missing_hpj_name"
        else:
            relation = classify_name_relation(spectrum.compound_name_raw, entry.name)
        relations.append(relation)
        if strict and relation in {"mismatch", "missing_hpj_name"}:
            raise ValueError(
                "Unsafe companion alignment at release row "
                f"{index}: HPJ={spectrum.compound_name_raw!r}, companion={entry.name!r}"
            )

    aligned: list[RawSpectrum] = []
    for spectrum, entry, relation in zip(spectra, entries, relations, strict=True):
        metadata = dict(spectrum.metadata)
        metadata.update(
            {
                "COMPOUND_LIST_NAME": entry.name,
                "COMPOUND_LIST_ROW": str(entry.row_number),
                "COMPOUND_LIST_NAME_RELATION": relation,
            }
        )
        if entry.new_in_source_version:
            metadata["COMPOUND_LIST_NEW_IN_VERSION"] = entry.new_in_source_version
        aligned.append(replace(spectrum, metadata=metadata))

    counts = Counter(relations)
    versions = {entry.new_in_source_version for entry in entries if entry.new_in_source_version}
    version = next(iter(versions)) if len(versions) == 1 else None
    report = CompanionAlignmentReport(
        number_entries=len(entries),
        relation_counts=dict(sorted(counts.items())),
        number_new_in_source_version=sum(
            entry.new_in_source_version is not None for entry in entries
        ),
        source_version_marker=version,
    )
    return aligned, report


def align_swgdrug_exact_mass_list(
    spectra: list[RawSpectrum],
    entries: list[ExactMassEntry],
    *,
    strict: bool = True,
) -> tuple[list[RawSpectrum], ExactMassAlignmentReport]:
    """Align the source formula/exact-mass workbook to spectra by release order."""
    if len(spectra) != len(entries):
        raise ValueError(
            "Spectrum/exact-mass row count mismatch: "
            f"{len(spectra)} spectra vs {len(entries)} exact-mass entries"
        )

    relations: list[str] = []
    formula_matches = 0
    aligned: list[RawSpectrum] = []
    for index, (spectrum, entry) in enumerate(zip(spectra, entries, strict=True), start=1):
        hpj_name = spectrum.compound_name_raw
        companion_name = spectrum.metadata.get("COMPOUND_LIST_NAME")
        if companion_name:
            relation = "exact" if companion_name == entry.name else "mismatch"
        elif hpj_name:
            relation = classify_name_relation(hpj_name, entry.name)
        else:
            relation = "missing_hpj_name"
        relations.append(relation)
        if strict and relation in {"mismatch", "missing_hpj_name"}:
            raise ValueError(
                "Unsafe exact-mass alignment at release row "
                f"{index}: HPJ={hpj_name!r}, workbook={entry.name!r}"
            )

        hpj_formula = spectrum.metadata.get("MOLFORM") or spectrum.metadata.get("FORMULA")
        formula_matches_this_row = hpj_formula == entry.formula
        if formula_matches_this_row:
            formula_matches += 1
        elif strict:
            raise ValueError(
                "Formula disagreement in exact-mass alignment at release row "
                f"{index}: HPJ={hpj_formula!r}, workbook={entry.formula!r}"
            )

        metadata = dict(spectrum.metadata)
        metadata.update(
            {
                "EXACT_MASS_LIST_ROW": str(entry.row_number),
                "EXACT_MASS_LIST_NAME": entry.name,
                "EXACT_MASS_LIST_NAME_RELATION": relation,
                "EXACT_MASS_LIST_FORMULA": entry.formula,
                "EXACT_MASS_LIST_FORMULA_RELATION": (
                    "exact" if formula_matches_this_row else "mismatch"
                ),
                "EXACT_MASS_LIST_MONOISOTOPIC_MASS": repr(entry.monoisotopic_mass),
            }
        )
        if entry.m_plus_proton_mass is not None:
            metadata["EXACT_MASS_LIST_M_PLUS_PROTON"] = repr(entry.m_plus_proton_mass)
        aligned.append(replace(spectrum, metadata=metadata))

    report = ExactMassAlignmentReport(
        number_entries=len(entries),
        name_relation_counts=dict(sorted(Counter(relations).items())),
        number_formula_matches=formula_matches,
        number_m_plus_proton_values=sum(entry.m_plus_proton_mass is not None for entry in entries),
        number_m_plus_proton_missing=sum(entry.m_plus_proton_mass is None for entry in entries),
        source_version_marker=None,
    )
    return aligned, report


@dataclass(frozen=True, slots=True)
class AgilentHeaderEntry:
    """Printable metadata strings embedded in one Agilent library header record."""

    record_number: int
    strings: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class AgilentHeaderAlignmentReport:
    """Alignment summary for SWGDRUG's Agilent-format HEADER metadata."""

    number_entries: int
    source_comment_relation_counts: dict[str, int]
    number_extended_source_comments: int
    number_short_names_omitted_from_printable_scan: int
    number_short_formulas_omitted_from_printable_scan: int
    archive_member: str
    header_sha256: str

    def to_dict(self) -> dict[str, object]:
        return {
            "number_entries": self.number_entries,
            "source_comment_relation_counts": dict(self.source_comment_relation_counts),
            "number_extended_source_comments": self.number_extended_source_comments,
            "number_short_names_omitted_from_printable_scan": (
                self.number_short_names_omitted_from_printable_scan
            ),
            "number_short_formulas_omitted_from_printable_scan": (
                self.number_short_formulas_omitted_from_printable_scan
            ),
            "archive_member": self.archive_member,
            "header_sha256": self.header_sha256,
        }


def _ascii_strings(data: bytes, *, minimum_length: int = 4) -> list[str]:
    """Return contiguous printable ASCII runs, similar to the Unix ``strings`` utility."""
    strings: list[str] = []
    current = bytearray()
    for value in data:
        if 32 <= value <= 126:
            current.append(value)
            continue
        if len(current) >= minimum_length:
            strings.append(current.decode("ascii"))
        current.clear()
    if len(current) >= minimum_length:
        strings.append(current.decode("ascii"))
    return strings


def read_swgdrug_agilent_header_archive(
    path: Path,
) -> tuple[list[AgilentHeaderEntry], str, str]:
    """Read printable metadata blocks from a SWGDRUG Agilent-format ZIP archive.

    This deliberately does not attempt to parse vendor spectral binary structures. It
    only reads the ASCII metadata embedded in the archive's ``HEADER`` member and
    requires downstream alignment against the already parsed JCAMP release.
    """
    import hashlib
    import zipfile

    with zipfile.ZipFile(path) as archive:
        candidates = [name for name in archive.namelist() if Path(name).name == "HEADER"]
        if len(candidates) != 1:
            raise ValueError(
                f"Expected exactly one Agilent HEADER member in {path}, found {candidates!r}"
            )
        member = candidates[0]
        data = archive.read(member)

    header_sha256 = hashlib.sha256(data).hexdigest()
    strings = _ascii_strings(data)
    markers = [index for index, value in enumerate(strings) if value.startswith("SWGDRUG 3.")]
    if not markers:
        raise ValueError(f"No SWGDRUG record markers found in Agilent HEADER member {member}")

    entries: list[AgilentHeaderEntry] = []
    for record_number, start in enumerate(markers, start=1):
        stop = markers[record_number] if record_number < len(markers) else len(strings)
        entries.append(
            AgilentHeaderEntry(
                record_number=record_number,
                strings=tuple(strings[start + 1 : stop]),
            )
        )
    return entries, member, header_sha256


def align_swgdrug_agilent_header(
    spectra: list[RawSpectrum],
    entries: list[AgilentHeaderEntry],
    *,
    archive_member: str,
    header_sha256: str,
    strict: bool = True,
) -> tuple[list[RawSpectrum], AgilentHeaderAlignmentReport]:
    """Use the Agilent header to restore full SWGDRUG source comments conservatively.

    The HPJ representation truncates many source comments. The Agilent-format header
    preserves longer printable comments. A comment is accepted only when it is either
    identical to the HPJ comment or has the HPJ comment as an exact prefix.
    """
    if len(spectra) != len(entries):
        raise ValueError(
            "Spectrum/Agilent-header record count mismatch: "
            f"{len(spectra)} spectra vs {len(entries)} header entries"
        )

    relations: list[str] = []
    short_names = 0
    short_formulas = 0
    aligned: list[RawSpectrum] = []

    for release_index, (spectrum, entry) in enumerate(zip(spectra, entries, strict=True), start=1):
        tokens = list(entry.strings)
        position = 0
        expected_names = {
            value
            for value in (
                spectrum.compound_name_raw,
                spectrum.metadata.get("COMPOUND_LIST_NAME"),
                spectrum.metadata.get("EXACT_MASS_LIST_NAME"),
            )
            if value
        }
        if position < len(tokens) and tokens[position] in expected_names:
            position += 1
        elif expected_names and all(len(name) < 4 for name in expected_names):
            short_names += 1
        elif strict:
            raise ValueError(
                "Agilent HEADER name alignment failed at release row "
                f"{release_index}: expected one of {sorted(expected_names)!r}, "
                f"found {tokens[position] if position < len(tokens) else None!r}"
            )

        formula = spectrum.metadata.get("MOLFORM") or spectrum.metadata.get("FORMULA")
        if formula and position < len(tokens) and tokens[position] == formula:
            position += 1
        elif formula and len(formula) < 4:
            short_formulas += 1
        elif formula and strict:
            raise ValueError(
                "Agilent HEADER formula alignment failed at release row "
                f"{release_index}: expected {formula!r}, "
                f"found {tokens[position] if position < len(tokens) else None!r}"
            )

        source_comment = tokens[position] if position < len(tokens) else None
        hpj_comment = spectrum.metadata.get("SAMPLE_DESCRIPTION")
        if source_comment is None:
            relation = "missing_agilent_comment"
        elif hpj_comment == source_comment:
            relation = "exact"
        elif hpj_comment and source_comment.startswith(hpj_comment):
            relation = "hpj_truncated"
        elif not hpj_comment:
            relation = "agilent_only"
        else:
            relation = "mismatch"

        if strict and relation in {"missing_agilent_comment", "mismatch"}:
            raise ValueError(
                "Agilent HEADER source-comment alignment failed at release row "
                f"{release_index}: HPJ={hpj_comment!r}, Agilent={source_comment!r}"
            )
        relations.append(relation)

        metadata = dict(spectrum.metadata)
        if source_comment:
            metadata["SOURCE_COMMENT_FULL"] = source_comment
            metadata["SOURCE_COMMENT_RELATION"] = relation
            metadata["SOURCE_COMMENT_SOURCE_FORMAT"] = "swgdrug_agilent_header"
        aligned.append(replace(spectrum, metadata=metadata))

    counts = Counter(relations)
    report = AgilentHeaderAlignmentReport(
        number_entries=len(entries),
        source_comment_relation_counts=dict(sorted(counts.items())),
        number_extended_source_comments=counts.get("hpj_truncated", 0),
        number_short_names_omitted_from_printable_scan=short_names,
        number_short_formulas_omitted_from_printable_scan=short_formulas,
        archive_member=archive_member,
        header_sha256=header_sha256,
    )
    return aligned, report
