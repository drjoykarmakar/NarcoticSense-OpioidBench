"""Conservative chemical identifier normalization and validation."""

from __future__ import annotations

import hashlib
import re
import unicodedata

_CAS_RE = re.compile(r"^(?P<first>\d{2,7})-(?P<middle>\d{2})-(?P<check>\d)$")
_INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
_FORMULA_RE = re.compile(r"^(?:[A-Z][a-z]?\d*)+$")
_FORMULA_TOKEN_RE = re.compile(r"([A-Z][a-z]?)(\d*)")
CAS_PLACEHOLDERS = frozenset({"0", "000000-00-0", "0-00-0"})


def normalize_compound_name(name: str | None) -> str | None:
    """Normalize Unicode and whitespace only; do not infer synonyms or merge names."""
    if name is None:
        return None
    normalized = unicodedata.normalize("NFKC", name)
    normalized = " ".join(normalized.split())
    return normalized or None


def cas_checksum_valid(value: str | None) -> bool:
    """Return whether a CAS Registry Number has a valid check digit."""
    if not value:
        return False
    compact = value.strip()
    if compact in CAS_PLACEHOLDERS:
        return False
    match = _CAS_RE.fullmatch(compact)
    if not match:
        return False
    digits = match.group("first") + match.group("middle")
    weighted = sum(int(digit) * weight for weight, digit in enumerate(reversed(digits), 1))
    return weighted % 10 == int(match.group("check"))


def normalize_cas(value: str | None) -> str | None:
    """Return canonical CAS formatting when valid; otherwise return ``None``."""
    if not value:
        return None
    compact = value.strip()
    if not cas_checksum_valid(compact):
        return None
    match = _CAS_RE.fullmatch(compact)
    assert match is not None
    first = match.group("first").lstrip("0") or "0"
    return f"{first}-{match.group('middle')}-{match.group('check')}"


def inchikey_valid(value: str | None) -> bool:
    if not value:
        return False
    return bool(_INCHIKEY_RE.fullmatch(value.strip().upper()))


def molecular_formula_plausible(value: str | None) -> bool:
    """Perform syntax-only formula validation; this is not a chemical plausibility proof."""
    if not value:
        return False
    return bool(_FORMULA_RE.fullmatch(value.strip()))


def formula_atom_counts(value: str | None) -> dict[str, int] | None:
    """Return simple formula atom counts, ignoring a terminal charge annotation.

    The source formulas in this project are flat Hill-style formulas rather than
    salts/mixtures with parentheses or dots. Unsupported syntax returns ``None``
    rather than being guessed.
    """
    if not value:
        return None
    text = re.sub(r"[+-]\d*$", "", value.strip())
    counts: dict[str, int] = {}
    position = 0
    for match in _FORMULA_TOKEN_RE.finditer(text):
        if match.start() != position:
            return None
        element = match.group(1)
        counts[element] = counts.get(element, 0) + int(match.group(2) or "1")
        position = match.end()
    if not counts or position != len(text):
        return None
    return counts


def structure_formula_relation(
    source_formula: str | None, structure_formula: str | None
) -> str | None:
    """Classify a structure formula as exact or an unlabeled deuterated parent.

    NIST molfiles frequently omit explicit deuterium labeling while preserving
    the same heavy-atom connectivity. That relation is retained explicitly so
    such structures can support scaffold/fingerprint work without pretending
    they encode the isotope labels in the source formula.
    """
    source_counts = formula_atom_counts(source_formula)
    structure_counts = formula_atom_counts(structure_formula)
    if source_counts is None or structure_counts is None:
        return None
    if source_counts == structure_counts:
        return "exact_formula"
    if source_counts.get("D", 0):
        parent_counts = dict(source_counts)
        parent_counts["H"] = parent_counts.get("H", 0) + parent_counts.pop("D")
        if parent_counts == structure_counts:
            return "isotopologue_parent"
    return None


def conservative_compound_id(*, cas: str | None, name: str | None) -> str:
    """Build a stable identity key without fuzzy name matching.

    Valid CAS is preferred.  Otherwise an exact normalized source name is hashed so that
    unrelated compounds are never merged merely because their names appear similar.
    """
    normalized_cas = normalize_cas(cas)
    if normalized_cas:
        return f"cas:{normalized_cas}"
    normalized_name = normalize_compound_name(name)
    if normalized_name:
        digest = hashlib.sha256(normalized_name.encode("utf-8")).hexdigest()[:20]
        return f"name:{digest}"
    raise ValueError("Cannot build compound identity without a valid CAS or compound name")
