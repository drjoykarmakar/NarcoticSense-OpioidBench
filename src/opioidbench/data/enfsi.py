"""Restricted ENFSI-DWG entry-list inspection and cross-source audit utilities.

The ENFSI Drugs Working Group publishes a data-usage agreement that restricts
redistribution.  This module therefore emits aggregate cross-check statistics only;
it does not copy ENFSI row-level names or other source data into public artifacts.
"""

from __future__ import annotations

import json
import re
import unicodedata
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True, slots=True)
class ENFSIEntry:
    """One row from the ENFSI-DWG MS library entry list."""

    source_index: int
    name: str
    formula: str | None
    molecular_weight: float


@dataclass(frozen=True, slots=True)
class ENFSICrosscheckReport:
    """Aggregate-only comparison between ENFSI metadata and a local OpioidBench build."""

    enfsi_entries: int
    enfsi_unique_names: int
    enfsi_unique_name_formula_pairs: int
    enfsi_missing_formula: int
    enfsi_exact_duplicate_rows: int
    high_confidence_matched_entries: int
    high_confidence_matched_unique_swg_compounds: int
    ambiguous_entries: int
    unmatched_entries: int

    def to_dict(self) -> dict[str, int]:
        return {
            "enfsi_entries": self.enfsi_entries,
            "enfsi_unique_names": self.enfsi_unique_names,
            "enfsi_unique_name_formula_pairs": self.enfsi_unique_name_formula_pairs,
            "enfsi_missing_formula": self.enfsi_missing_formula,
            "enfsi_exact_duplicate_rows": self.enfsi_exact_duplicate_rows,
            "high_confidence_matched_entries": self.high_confidence_matched_entries,
            "high_confidence_matched_unique_swg_compounds": (
                self.high_confidence_matched_unique_swg_compounds
            ),
            "ambiguous_entries": self.ambiguous_entries,
            "unmatched_entries": self.unmatched_entries,
        }


def _load_openpyxl():
    try:
        from openpyxl import load_workbook
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise RuntimeError(
            "Reading ENFSI XLSX entry lists requires the optional 'excel' dependency. "
            "Install opioidbench[excel]."
        ) from exc
    return load_workbook


def read_enfsi_entry_list(path: Path) -> list[ENFSIEntry]:
    """Read an ENFSI-DWG MS entry-list workbook without re-exporting row-level data."""
    load_workbook = _load_openpyxl()
    workbook = load_workbook(path, read_only=True, data_only=True)
    if not workbook.sheetnames:
        raise ValueError(f"Workbook {path} contains no worksheets")
    sheet = workbook[workbook.sheetnames[0]]
    header = [sheet.cell(1, column).value for column in range(1, 5)]
    expected = ["#", "Name", "Formula", "MW"]
    if header != expected:
        raise ValueError(f"Unexpected ENFSI entry-list header: {header!r}")

    entries: list[ENFSIEntry] = []
    for row_number, row in enumerate(
        sheet.iter_rows(min_row=2, max_col=4, values_only=True), start=2
    ):
        source_index, name, formula, molecular_weight = row
        if all(value is None for value in row):
            continue
        if source_index is None or name is None or molecular_weight is None:
            raise ValueError(f"Incomplete ENFSI entry-list row {row_number}")
        try:
            parsed_index = int(source_index)
            parsed_weight = float(molecular_weight)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Invalid ENFSI numeric value at row {row_number}") from exc
        if parsed_weight <= 0:
            raise ValueError(f"Non-positive ENFSI molecular weight at row {row_number}")
        entries.append(
            ENFSIEntry(
                source_index=parsed_index,
                name=str(name).strip(),
                formula=(
                    str(formula).strip() if formula is not None and str(formula).strip() else None
                ),
                molecular_weight=parsed_weight,
            )
        )
    return entries


def _normalize_name(value: str) -> str:
    ascii_value = unicodedata.normalize("NFKD", value).encode("ascii", "ignore").decode()
    return re.sub(r"[^a-z0-9]+", "", ascii_value.casefold())


def _conservative_aliases(name: str) -> tuple[str, ...]:
    """Extract explicit aliases from source punctuation without spelling corrections."""
    aliases: list[str] = [name]
    aliases.extend(part.strip() for part in name.split(";") if part.strip())
    aliases.extend(
        match.group(1).strip()
        for match in re.finditer(r"\(([^()]*)\)", name)
        if match.group(1).strip()
    )
    outside = re.sub(r"\([^()]*\)", " ", name).strip()
    if outside:
        aliases.append(outside)
    return tuple(dict.fromkeys(aliases))


def _load_normalized_records(path: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                payload = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(
                    f"Invalid normalized-record JSON at line {line_number} in {path}"
                ) from exc
            if not isinstance(payload, dict):
                raise ValueError(f"Normalized record at line {line_number} must be an object")
            records.append(payload)
    return records


def audit_enfsi_entry_list(
    entry_list_path: Path,
    normalized_records_path: Path,
) -> ENFSICrosscheckReport:
    """Cross-check ENFSI metadata against a local SWGDRUG-derived normalized build.

    A high-confidence match requires an explicit ENFSI alias, molecular formula, and
    rounded nominal molecular weight to resolve to exactly one SWGDRUG spectrum record.
    The report intentionally contains aggregate counts only.
    """
    entries = read_enfsi_entry_list(entry_list_path)
    records = _load_normalized_records(normalized_records_path)

    index: dict[tuple[str, str | None, int | None], set[str]] = defaultdict(set)
    for record in records:
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
                    (
                        _normalize_name(name),
                        record.get("molecular_formula"),
                        rounded_weight,
                    )
                ].add(str(record.get("compound_id") or record.get("spectrum_id") or ""))

    matched_entries = 0
    matched_compound_ids: set[str] = set()
    ambiguous_entries = 0
    unmatched_entries = 0
    for entry in entries:
        candidates: set[int] = set()
        rounded_weight = round(entry.molecular_weight)
        for alias in _conservative_aliases(entry.name):
            candidates.update(
                index.get((_normalize_name(alias), entry.formula, rounded_weight), set())
            )
        if len(candidates) == 1:
            matched_entries += 1
            matched_compound_ids.update(candidates)
        elif len(candidates) > 1:
            ambiguous_entries += 1
        else:
            unmatched_entries += 1

    exact_rows = Counter((entry.name, entry.formula, entry.molecular_weight) for entry in entries)
    return ENFSICrosscheckReport(
        enfsi_entries=len(entries),
        enfsi_unique_names=len({entry.name for entry in entries}),
        enfsi_unique_name_formula_pairs=len({(entry.name, entry.formula) for entry in entries}),
        enfsi_missing_formula=sum(entry.formula is None for entry in entries),
        enfsi_exact_duplicate_rows=sum(count - 1 for count in exact_rows.values() if count > 1),
        high_confidence_matched_entries=matched_entries,
        high_confidence_matched_unique_swg_compounds=len(matched_compound_ids),
        ambiguous_entries=ambiguous_entries,
        unmatched_entries=unmatched_entries,
    )
