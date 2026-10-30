"""Evidence-backed taxonomy curation without keyword or similarity inference."""

from __future__ import annotations

import csv
import hashlib
from dataclasses import dataclass
from pathlib import Path

from opioidbench.chemistry.exclusions import (
    analytical_variant_reason,
    load_taxonomy_exclusions,
    match_taxonomy_exclusion,
)
from opioidbench.chemistry.identifiers import normalize_cas, normalize_compound_name
from opioidbench.chemistry.taxonomy import REQUIRED_COLUMNS

EVIDENCE_COLUMNS = (
    "evidence_id",
    "preferred_name",
    "aliases",
    "cas",
    "opioid_status",
    "opioid_family",
    "subfamily",
    "evidence_source",
    "evidence_url",
    "evidence_date",
    "identity_source",
    "identity_url",
    "notes",
)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _name_key(value: str | None) -> str | None:
    normalized = normalize_compound_name(value)
    return normalized.casefold() if normalized else None


def _split_aliases(value: str | None) -> tuple[str, ...]:
    if not value:
        return ()
    return tuple(part.strip() for part in value.split("|") if part.strip())


@dataclass(frozen=True, slots=True)
class OpioidEvidence:
    evidence_id: str
    preferred_name: str
    aliases: tuple[str, ...]
    cas: str | None
    opioid_status: str
    opioid_family: str
    subfamily: str
    evidence_source: str
    evidence_url: str
    evidence_date: str
    identity_source: str
    identity_url: str
    notes: str

    @property
    def name_keys(self) -> frozenset[str]:
        values = (self.preferred_name, *self.aliases)
        keys = {_name_key(value) for value in values}
        return frozenset(key for key in keys if key)


def load_opioid_evidence(path: Path) -> list[OpioidEvidence]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = [column for column in EVIDENCE_COLUMNS if column not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"Opioid evidence file is missing columns: {', '.join(missing)}")
        rows: list[OpioidEvidence] = []
        seen_ids: set[str] = set()
        for raw in reader:
            evidence_id = (raw.get("evidence_id") or "").strip()
            if not evidence_id:
                raise ValueError("Opioid evidence rows require evidence_id")
            if evidence_id in seen_ids:
                raise ValueError(f"Duplicate opioid evidence_id: {evidence_id}")
            seen_ids.add(evidence_id)
            status = (raw.get("opioid_status") or "").strip().lower()
            if status not in {"opioid", "non-opioid"}:
                raise ValueError(
                    f"Evidence row {evidence_id!r} has unsupported opioid_status {status!r}"
                )
            preferred_name = (raw.get("preferred_name") or "").strip()
            if not preferred_name:
                raise ValueError(f"Evidence row {evidence_id!r} requires preferred_name")
            cas_raw = (raw.get("cas") or "").strip()
            cas = normalize_cas(cas_raw) if cas_raw else None
            if cas_raw and not cas:
                raise ValueError(f"Evidence row {evidence_id!r} has invalid CAS {cas_raw!r}")
            rows.append(
                OpioidEvidence(
                    evidence_id=evidence_id,
                    preferred_name=preferred_name,
                    aliases=_split_aliases(raw.get("aliases")),
                    cas=cas,
                    opioid_status=status,
                    opioid_family=(raw.get("opioid_family") or "").strip(),
                    subfamily=(raw.get("subfamily") or "").strip(),
                    evidence_source=(raw.get("evidence_source") or "").strip(),
                    evidence_url=(raw.get("evidence_url") or "").strip(),
                    evidence_date=(raw.get("evidence_date") or "").strip(),
                    identity_source=(raw.get("identity_source") or "").strip(),
                    identity_url=(raw.get("identity_url") or "").strip(),
                    notes=(raw.get("notes") or "").strip(),
                )
            )
    return rows


def _compound_name_keys(row: dict[str, str]) -> set[str]:
    values = [row.get("preferred_name", ""), *_split_aliases(row.get("aliases"))]
    keys = {_name_key(value) for value in values}
    return {key for key in keys if key}


def _match_evidence(
    row: dict[str, str], evidence: list[OpioidEvidence]
) -> list[tuple[OpioidEvidence, str]]:
    cas = normalize_cas(row.get("cas"))
    names = _compound_name_keys(row)
    matches: list[tuple[OpioidEvidence, str]] = []
    for item in evidence:
        if cas and item.cas and cas == item.cas:
            matches.append((item, "cas"))
            continue
        if names & set(item.name_keys):
            matches.append((item, "name"))
    return matches


def curate_taxonomy_from_evidence(
    compounds_path: Path,
    evidence_path: Path,
    *,
    taxonomy_output: Path,
    review_output: Path,
    exclusions_path: Path | None = None,
) -> dict[str, int | str]:
    """Write exact curated labels; send unsupported/conflicting compounds to review."""
    evidence = load_opioid_evidence(evidence_path)
    exclusions = load_taxonomy_exclusions(exclusions_path) if exclusions_path else []
    with compounds_path.open(newline="", encoding="utf-8") as handle:
        compounds = list(csv.DictReader(handle))

    taxonomy_rows: list[dict[str, str]] = []
    review_rows: list[dict[str, str]] = []
    counts: dict[str, int | str] = {
        "input_compounds": len(compounds),
        "input_compounds_sha256": _sha256_file(compounds_path),
        "evidence_records": len(evidence),
        "evidence_sha256": _sha256_file(evidence_path),
        "exclusion_records": len(exclusions),
        "exclusions_sha256": _sha256_file(exclusions_path) if exclusions_path else "",
        "curated_opioid": 0,
        "curated_nonopioid": 0,
        "uncertain_conflict": 0,
        "excluded_target": 0,
        "excluded_analytical_variant": 0,
        "pending": 0,
        "matched_by_cas": 0,
        "matched_by_name": 0,
    }

    for row in compounds:
        preferred_name = (row.get("preferred_name") or "").strip()
        compound_id = (row.get("compound_id") or "").strip()
        exclusion = match_taxonomy_exclusion(row, exclusions) if exclusions else None
        if exclusion:
            counts["excluded_target"] += 1
            taxonomy_rows.append(
                {
                    "compound_id": compound_id,
                    "preferred_name": preferred_name,
                    "opioid_status": "",
                    "opioid_family": "benchmark exclusion",
                    "subfamily": exclusion.reason,
                    "evidence_source": exclusion.evidence_source,
                    "curation_status": "excluded",
                    "notes": f"Target exclusion; exclusion_id={exclusion.exclusion_id}",
                }
            )
            continue
        analytical_reason = analytical_variant_reason(row)
        if analytical_reason:
            counts["excluded_analytical_variant"] += 1
            taxonomy_rows.append(
                {
                    "compound_id": compound_id,
                    "preferred_name": preferred_name,
                    "opioid_status": "",
                    "opioid_family": "benchmark exclusion",
                    "subfamily": analytical_reason,
                    "evidence_source": "OpioidBench analytical-variant policy",
                    "curation_status": "excluded",
                    "notes": "Measurement standard/preparation excluded from biological labels.",
                }
            )
            continue
        matches = _match_evidence(row, evidence)
        if not matches:
            counts["pending"] += 1
            review_rows.append(
                {
                    "compound_id": compound_id,
                    "preferred_name": preferred_name,
                    "aliases": row.get("aliases", ""),
                    "cas": row.get("cas", ""),
                    "molecular_formula": row.get("molecular_formula", ""),
                    "inchikey": row.get("inchikey", ""),
                    "structure_available": "yes" if row.get("canonical_smiles") else "no",
                    "review_reason": "no_exact_authoritative_evidence_match",
                    "recommended_action": (
                        "Review WHO ECDD/ATC or another traceable pharmacology source; "
                        "do not infer from the compound name alone."
                    ),
                }
            )
            continue

        statuses = {item.opioid_status for item, _ in matches}
        methods = {method for _, method in matches}
        if len(statuses) != 1:
            counts["uncertain_conflict"] += 1
            review_rows.append(
                {
                    "compound_id": compound_id,
                    "preferred_name": preferred_name,
                    "aliases": row.get("aliases", ""),
                    "cas": row.get("cas", ""),
                    "molecular_formula": row.get("molecular_formula", ""),
                    "inchikey": row.get("inchikey", ""),
                    "structure_available": "yes" if row.get("canonical_smiles") else "no",
                    "review_reason": "conflicting_authoritative_evidence",
                    "recommended_action": "Resolve the evidence conflict manually.",
                }
            )
            continue

        status = next(iter(statuses))
        families = {item.opioid_family for item, _ in matches if item.opioid_family}
        subfamilies = {item.subfamily for item, _ in matches if item.subfamily}
        sources = {item.evidence_source for item, _ in matches if item.evidence_source}
        sources.update(item.identity_source for item, _ in matches if item.identity_source)
        sources = sorted(sources)
        evidence_ids = sorted({item.evidence_id for item, _ in matches})
        family = next(iter(families)) if len(families) == 1 else " | ".join(sorted(families))
        subfamily = (
            next(iter(subfamilies)) if len(subfamilies) == 1 else " | ".join(sorted(subfamilies))
        )
        taxonomy_rows.append(
            {
                "compound_id": compound_id,
                "preferred_name": preferred_name,
                "opioid_status": status,
                "opioid_family": family,
                "subfamily": subfamily,
                "evidence_source": " | ".join(sources),
                "curation_status": "curated",
                "notes": (
                    f"Exact evidence match ({'+'.join(sorted(methods))}); evidence_id="
                    + "|".join(evidence_ids)
                ),
            }
        )
        counts[f"curated_{'nonopioid' if status == 'non-opioid' else status}"] += 1
        if "cas" in methods:
            counts["matched_by_cas"] += 1
        if "name" in methods:
            counts["matched_by_name"] += 1

    taxonomy_output.parent.mkdir(parents=True, exist_ok=True)
    with taxonomy_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(REQUIRED_COLUMNS), lineterminator="\n")
        writer.writeheader()
        writer.writerows(taxonomy_rows)

    review_output.parent.mkdir(parents=True, exist_ok=True)
    review_fields = [
        "compound_id",
        "preferred_name",
        "aliases",
        "cas",
        "molecular_formula",
        "inchikey",
        "structure_available",
        "review_reason",
        "recommended_action",
    ]
    with review_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=review_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(review_rows)

    return counts
