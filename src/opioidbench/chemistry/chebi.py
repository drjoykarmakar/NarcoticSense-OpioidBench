"""Streaming ChEBI bulk-SDF parsing and conservative compound matching."""

from __future__ import annotations

import csv
import gzip
import re
from collections import defaultdict
from collections.abc import Iterable, Iterator
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, TextIO

from opioidbench.chemistry.identifiers import (
    normalize_cas,
    normalize_compound_name,
    structure_formula_relation,
)
from opioidbench.chemistry.structures import calculate_structure_features

_CAS_TOKEN_RE = re.compile(r"\b\d{2,7}-\d{2}-\d\b")
_PROPERTY_RE = re.compile(r"^>\s*<([^>]+)>")


def _property_key(value: str) -> str:
    return re.sub(r"[^A-Z0-9]+", " ", value.upper()).strip()


def _name_key(value: str | None) -> str | None:
    normalized = normalize_compound_name(value)
    return normalized.casefold() if normalized else None


def _first_property(properties: dict[str, tuple[str, ...]], *keys: str) -> str | None:
    for key in keys:
        values = properties.get(_property_key(key))
        if values:
            for value in values:
                if value.strip():
                    return value.strip()
    return None


def _split_aliases(value: str | None) -> tuple[str, ...]:
    if not value:
        return ()
    return tuple(part.strip() for part in value.split("|") if part.strip())


@dataclass(frozen=True, slots=True)
class ChEBIEntry:
    chebi_id: str
    name: str | None
    synonyms: tuple[str, ...]
    cas_numbers: tuple[str, ...]
    formula: str | None
    monoisotopic_mass: float | None
    smiles: str | None
    inchi: str | None
    inchikey: str | None
    stars: int | None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class ChEBIMatch:
    compound_id: str
    status: str
    chebi_id: str | None
    match_method: str | None
    candidate_count: int
    canonical_smiles: str | None
    inchikey: str | None
    formula_relation: str | None
    error: str | None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _open_text(path: Path) -> TextIO:
    if path.suffix.lower() == ".gz":
        return gzip.open(path, mode="rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def _parse_sdf_properties(lines: list[str]) -> dict[str, tuple[str, ...]]:
    properties: dict[str, list[str]] = defaultdict(list)
    idx = 0
    while idx < len(lines):
        match = _PROPERTY_RE.match(lines[idx].rstrip("\n"))
        if not match:
            idx += 1
            continue
        key = _property_key(match.group(1))
        idx += 1
        values: list[str] = []
        while idx < len(lines):
            line = lines[idx].rstrip("\r\n")
            if not line:
                break
            if _PROPERTY_RE.match(line):
                break
            values.append(line)
            idx += 1
        properties[key].extend(values)
        while idx < len(lines) and not lines[idx].strip():
            idx += 1
    return {key: tuple(values) for key, values in properties.items()}


def _float_or_none(value: str | None) -> float | None:
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _int_or_none(value: str | None) -> int | None:
    if not value:
        return None
    try:
        return int(value)
    except ValueError:
        return None


def _entry_from_properties(properties: dict[str, tuple[str, ...]]) -> ChEBIEntry | None:
    chebi_id = _first_property(properties, "ChEBI ID", "CHEBI ID")
    if not chebi_id:
        return None
    name = _first_property(properties, "ChEBI NAME", "CHEBI NAME", "NAME")

    synonym_values: list[str] = []
    for key, values in properties.items():
        if "SYNONYM" in key:
            synonym_values.extend(value.strip() for value in values if value.strip())

    cas_numbers: set[str] = set()
    for key, values in properties.items():
        if key in {"CAS", "CAS REGISTRY NUMBER", "CAS REGISTRY NUMBERS"} or key.startswith("CAS "):
            for value in values:
                for token in _CAS_TOKEN_RE.findall(value):
                    normalized = normalize_cas(token)
                    if normalized:
                        cas_numbers.add(normalized)

    formula = _first_property(properties, "FORMULA", "MOLECULAR FORMULA")
    monoisotopic_mass = _float_or_none(
        _first_property(properties, "MONOISOTOPIC MASS", "MONOISOTOPICMASS")
    )
    smiles = _first_property(properties, "SMILES")
    inchi = _first_property(properties, "INCHI", "INCHI STRING")
    inchikey = _first_property(properties, "INCHIKEY", "INCHI KEY", "INCHI KEY STRING")
    stars = _int_or_none(_first_property(properties, "STAR", "STARS"))
    return ChEBIEntry(
        chebi_id=chebi_id,
        name=name,
        synonyms=tuple(dict.fromkeys(synonym_values)),
        cas_numbers=tuple(sorted(cas_numbers)),
        formula=formula,
        monoisotopic_mass=monoisotopic_mass,
        smiles=smiles,
        inchi=inchi,
        inchikey=inchikey.upper() if inchikey else None,
        stars=stars,
    )


def iter_chebi_sdf(path: Path) -> Iterator[ChEBIEntry]:
    """Stream ChEBI FULL/LITE SDF records without loading the bulk file into memory."""
    with _open_text(path) as handle:
        record: list[str] = []
        for line in handle:
            if line.rstrip("\r\n") == "$$$$":
                properties = _parse_sdf_properties(record)
                entry = _entry_from_properties(properties)
                if entry is not None:
                    yield entry
                record = []
            else:
                record.append(line)
        if record:
            properties = _parse_sdf_properties(record)
            entry = _entry_from_properties(properties)
            if entry is not None:
                yield entry


def _compound_name_keys(row: dict[str, str]) -> set[str]:
    values = [row.get("preferred_name", ""), *_split_aliases(row.get("aliases"))]
    keys = {_name_key(value) for value in values}
    return {key for key in keys if key}


def _entry_name_keys(entry: ChEBIEntry) -> set[str]:
    keys = {_name_key(entry.name), *(_name_key(value) for value in entry.synonyms)}
    return {key for key in keys if key}


def _index_candidate_rows(rows: list[dict[str, str]]):
    by_cas: dict[str, set[int]] = defaultdict(set)
    by_inchikey: dict[str, set[int]] = defaultdict(set)
    by_name: dict[str, set[int]] = defaultdict(set)
    for idx, row in enumerate(rows):
        cas = normalize_cas(row.get("cas"))
        if cas:
            by_cas[cas].add(idx)
        inchikey = (row.get("inchikey") or "").strip().upper()
        if inchikey:
            by_inchikey[inchikey].add(idx)
        for key in _compound_name_keys(row):
            by_name[key].add(idx)
    return by_cas, by_inchikey, by_name


def collect_relevant_chebi_entries(
    path: Path, rows: list[dict[str, str]]
) -> tuple[dict[int, list[tuple[ChEBIEntry, tuple[str, ...]]]], int]:
    """Scan a ChEBI SDF once and retain only entries touching exact target keys."""
    by_cas, by_inchikey, by_name = _index_candidate_rows(rows)
    candidates: dict[int, dict[str, tuple[ChEBIEntry, set[str]]]] = defaultdict(dict)
    scanned = 0
    for entry in iter_chebi_sdf(path):
        scanned += 1
        row_methods: dict[int, set[str]] = defaultdict(set)
        if entry.inchikey and entry.inchikey in by_inchikey:
            for idx in by_inchikey[entry.inchikey]:
                row_methods[idx].add("inchikey")
        for cas in entry.cas_numbers:
            for idx in by_cas.get(cas, ()):
                row_methods[idx].add("cas")
        for name_key in _entry_name_keys(entry):
            for idx in by_name.get(name_key, ()):
                row_methods[idx].add("name")
        for idx, methods in row_methods.items():
            previous = candidates[idx].get(entry.chebi_id)
            if previous is None:
                candidates[idx][entry.chebi_id] = (entry, set(methods))
            else:
                previous[1].update(methods)
    output: dict[int, list[tuple[ChEBIEntry, tuple[str, ...]]]] = {}
    for idx, per_id in candidates.items():
        output[idx] = [
            (entry, tuple(sorted(methods)))
            for entry, methods in sorted(per_id.values(), key=lambda item: item[0].chebi_id)
        ]
    return output, scanned


def _match_rank(methods: Iterable[str]) -> int:
    method_set = set(methods)
    if "inchikey" in method_set:
        return 3
    if "cas" in method_set:
        return 2
    if "name" in method_set:
        return 1
    return 0


def _is_true(value: str | None) -> bool:
    return (value or "").strip().casefold() in {"1", "true", "yes"}


def _validate_chebi_candidate(
    row: dict[str, str], entry: ChEBIEntry, methods: tuple[str, ...]
) -> tuple[bool, str | None, str | None, str | None, str | None]:
    """Return accepted, canonical SMILES, InChIKey, formula relation, error."""
    source_conflict = any(
        _is_true(row.get(field))
        for field in (
            "name_conflict",
            "formula_conflict",
            "exact_mass_conflict",
            "structure_conflict",
        )
    )
    if source_conflict and "inchikey" not in methods:
        return False, None, None, None, "source_identity_conflict"

    source_formula = (row.get("molecular_formula") or "").strip() or None
    formula_relation = None
    if source_formula and entry.formula:
        formula_relation = structure_formula_relation(source_formula, entry.formula)
        if formula_relation is None:
            return False, None, None, None, "source_formula_mismatch"
    elif "name" in methods and not ({"cas", "inchikey"} & set(methods)):
        return False, None, None, None, "name_match_missing_formula_evidence"

    canonical_smiles = None
    external_inchikey = entry.inchikey
    if entry.smiles:
        try:
            features = calculate_structure_features(entry.smiles)
        except ValueError:
            return False, None, None, formula_relation, "invalid_chebi_smiles"
        canonical_smiles = features.canonical_smiles
        if external_inchikey and features.inchikey != external_inchikey:
            return False, None, None, formula_relation, "chebi_smiles_inchikey_mismatch"
        external_inchikey = features.inchikey

    existing_inchikey = (row.get("inchikey") or "").strip().upper() or None
    if existing_inchikey and external_inchikey and existing_inchikey != external_inchikey:
        return False, None, None, formula_relation, "existing_structure_conflict"

    source_exact_mass = (row.get("source_monoisotopic_mass") or "").strip()
    if source_exact_mass and entry.monoisotopic_mass is not None:
        try:
            delta = abs(float(source_exact_mass) - entry.monoisotopic_mass)
        except ValueError:
            return False, None, None, formula_relation, "invalid_source_exact_mass"
        if delta > 0.0001 and formula_relation == "exact_formula":
            return False, None, None, formula_relation, "source_exact_mass_mismatch"

    return True, canonical_smiles, external_inchikey, formula_relation, None


def select_chebi_match(
    row: dict[str, str], candidates: list[tuple[ChEBIEntry, tuple[str, ...]]]
) -> ChEBIMatch:
    """Select only a unique, consistency-checked highest-confidence ChEBI candidate."""
    compound_id = row.get("compound_id", "")
    accepted: list[tuple[int, ChEBIEntry, tuple[str, ...], str | None, str | None, str | None]] = []
    errors: list[str] = []
    for entry, methods in candidates:
        ok, smiles, inchikey, formula_relation, error = _validate_chebi_candidate(
            row, entry, methods
        )
        if ok:
            accepted.append(
                (_match_rank(methods), entry, methods, smiles, inchikey, formula_relation)
            )
        elif error:
            errors.append(f"{entry.chebi_id}:{error}")

    if not accepted:
        status = "conflict" if candidates else "unmatched"
        return ChEBIMatch(
            compound_id=compound_id,
            status=status,
            chebi_id=None,
            match_method=None,
            candidate_count=len(candidates),
            canonical_smiles=None,
            inchikey=None,
            formula_relation=None,
            error=" | ".join(errors) if errors else None,
        )

    top_rank = max(item[0] for item in accepted)
    top = [item for item in accepted if item[0] == top_rank]
    if len(top) != 1:
        return ChEBIMatch(
            compound_id=compound_id,
            status="ambiguous",
            chebi_id=None,
            match_method=None,
            candidate_count=len(top),
            canonical_smiles=None,
            inchikey=None,
            formula_relation=None,
            error="multiple_consistent_chebi_candidates",
        )

    _, entry, methods, smiles, inchikey, formula_relation = top[0]
    method = "+".join(name for name in ("inchikey", "cas", "name") if name in set(methods))
    return ChEBIMatch(
        compound_id=compound_id,
        status="resolved",
        chebi_id=entry.chebi_id,
        match_method=method,
        candidate_count=len(candidates),
        canonical_smiles=smiles,
        inchikey=inchikey,
        formula_relation=formula_relation,
        error=None,
    )


def resolve_chebi_compound_table(
    input_path: Path,
    *,
    chebi_sdf: Path,
    output_path: Path | None = None,
) -> tuple[list[ChEBIMatch], dict[str, int]]:
    """Resolve a compound inventory against a local ChEBI FULL/LITE SDF dump."""
    with input_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    candidates, scanned = collect_relevant_chebi_entries(chebi_sdf, rows)
    matches = [select_chebi_match(row, candidates.get(idx, [])) for idx, row in enumerate(rows)]
    counts = {
        "input_compounds": len(rows),
        "chebi_entries_scanned": scanned,
        "resolved": sum(item.status == "resolved" for item in matches),
        "ambiguous": sum(item.status == "ambiguous" for item in matches),
        "conflict": sum(item.status == "conflict" for item in matches),
        "unmatched": sum(item.status == "unmatched" for item in matches),
        "structures_available_from_chebi": sum(
            item.status == "resolved" and item.canonical_smiles is not None for item in matches
        ),
    }
    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fieldnames = list(ChEBIMatch.__dataclass_fields__)
        with output_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(item.to_dict() for item in matches)
    return matches, counts
