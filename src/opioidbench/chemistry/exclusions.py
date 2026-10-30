"""Exact benchmark-target exclusions and analytical-variant detection."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from opioidbench.chemistry.identifiers import normalize_cas, normalize_compound_name

EXCLUSION_COLUMNS = (
    "exclusion_id",
    "preferred_name",
    "aliases",
    "cas",
    "reason",
    "evidence_source",
    "evidence_url",
    "evidence_date",
    "notes",
)


def _name_key(value: str | None) -> str | None:
    normalized = normalize_compound_name(value)
    return normalized.casefold() if normalized else None


def _split_aliases(value: str | None) -> tuple[str, ...]:
    if not value:
        return ()
    return tuple(part.strip() for part in value.split("|") if part.strip())


@dataclass(frozen=True, slots=True)
class TaxonomyExclusion:
    exclusion_id: str
    preferred_name: str
    aliases: tuple[str, ...]
    cas: str | None
    reason: str
    evidence_source: str
    evidence_url: str
    evidence_date: str
    notes: str

    @property
    def name_keys(self) -> frozenset[str]:
        values = (self.preferred_name, *self.aliases)
        keys = {_name_key(value) for value in values}
        return frozenset(key for key in keys if key)


def load_taxonomy_exclusions(path: Path) -> list[TaxonomyExclusion]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = [
            column for column in EXCLUSION_COLUMNS if column not in (reader.fieldnames or [])
        ]
        if missing:
            raise ValueError(f"Taxonomy exclusion file is missing columns: {', '.join(missing)}")
        rows: list[TaxonomyExclusion] = []
        seen_ids: set[str] = set()
        for raw in reader:
            exclusion_id = (raw.get("exclusion_id") or "").strip()
            if not exclusion_id:
                raise ValueError("Taxonomy exclusion rows require exclusion_id")
            if exclusion_id in seen_ids:
                raise ValueError(f"Duplicate exclusion_id: {exclusion_id}")
            seen_ids.add(exclusion_id)
            preferred_name = (raw.get("preferred_name") or "").strip()
            if not preferred_name:
                raise ValueError(f"Exclusion row {exclusion_id!r} requires preferred_name")
            cas_raw = (raw.get("cas") or "").strip()
            cas = normalize_cas(cas_raw) if cas_raw else None
            if cas_raw and not cas:
                raise ValueError(f"Exclusion row {exclusion_id!r} has invalid CAS {cas_raw!r}")
            reason = (raw.get("reason") or "").strip()
            if not reason:
                raise ValueError(f"Exclusion row {exclusion_id!r} requires reason")
            rows.append(
                TaxonomyExclusion(
                    exclusion_id=exclusion_id,
                    preferred_name=preferred_name,
                    aliases=_split_aliases(raw.get("aliases")),
                    cas=cas,
                    reason=reason,
                    evidence_source=(raw.get("evidence_source") or "").strip(),
                    evidence_url=(raw.get("evidence_url") or "").strip(),
                    evidence_date=(raw.get("evidence_date") or "").strip(),
                    notes=(raw.get("notes") or "").strip(),
                )
            )
    return rows


def match_taxonomy_exclusion(
    row: dict[str, str], exclusions: list[TaxonomyExclusion]
) -> TaxonomyExclusion | None:
    cas = normalize_cas(row.get("cas"))
    values = [row.get("preferred_name", ""), *_split_aliases(row.get("aliases"))]
    names = {_name_key(value) for value in values}
    names.discard(None)
    matches: list[TaxonomyExclusion] = []
    for item in exclusions:
        if cas and item.cas and cas == item.cas:
            matches.append(item)
            continue
        if names & set(item.name_keys):
            matches.append(item)
    if len(matches) > 1:
        ids = ", ".join(sorted(item.exclusion_id for item in matches))
        raise ValueError(f"Compound matches multiple taxonomy exclusions: {ids}")
    return matches[0] if matches else None


def analytical_variant_reason(row: dict[str, str]) -> str | None:
    """Return a deterministic benchmark-exclusion reason for analytical variants.

    These checks identify measurement standards/preparations, not biological drug classes.
    """
    formula = (row.get("molecular_formula") or "").strip()
    smiles = (row.get("canonical_smiles") or "").strip()
    name = normalize_compound_name(row.get("preferred_name")) or ""
    # Deuterium is represented as D in source formulae and/or [2H] in canonical SMILES.
    import re

    if "[2H]" in smiles or re.search(r"(?:^|(?<=[0-9]))D[0-9]*(?=[A-Z]|$)", formula):
        return "isotope_labeled_analytical_standard"
    if "tms" in name.casefold():
        return "derivatized_analytical_standard"
    return None
