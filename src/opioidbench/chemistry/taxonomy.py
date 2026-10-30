"""Human-reviewed opioid taxonomy schema and loading utilities."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

REQUIRED_COLUMNS = (
    "compound_id",
    "preferred_name",
    "opioid_status",
    "opioid_family",
    "subfamily",
    "evidence_source",
    "curation_status",
    "notes",
)
ALLOWED_CURATION_STATUS = frozenset({"pending", "curated", "uncertain", "excluded"})


@dataclass(frozen=True, slots=True)
class TaxonomyEntry:
    compound_id: str
    preferred_name: str
    opioid_status: str
    opioid_family: str
    subfamily: str
    evidence_source: str
    curation_status: str
    notes: str


def load_taxonomy(path: Path) -> dict[str, TaxonomyEntry]:
    """Load a curated taxonomy table without keyword-based opioid inference."""
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != REQUIRED_COLUMNS:
            raise ValueError(
                f"Taxonomy columns must be exactly {REQUIRED_COLUMNS}; got {reader.fieldnames}"
            )
        entries: dict[str, TaxonomyEntry] = {}
        for row_number, row in enumerate(reader, start=2):
            entry = TaxonomyEntry(**{column: row[column].strip() for column in REQUIRED_COLUMNS})
            if not entry.compound_id:
                raise ValueError(f"Missing compound_id at taxonomy row {row_number}")
            if entry.curation_status not in ALLOWED_CURATION_STATUS:
                raise ValueError(
                    f"Invalid curation_status {entry.curation_status!r} at row {row_number}"
                )
            if entry.compound_id in entries:
                raise ValueError(f"Duplicate taxonomy compound_id {entry.compound_id!r}")
            entries[entry.compound_id] = entry
    return entries
