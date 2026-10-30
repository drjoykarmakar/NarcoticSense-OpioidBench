"""Conservative ChEBI-first, PubChem-second chemical identity resolution."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any

from opioidbench.chemistry.chebi import ChEBIMatch, resolve_chebi_compound_table
from opioidbench.chemistry.enrichment import EnrichmentRecord, MetadataCache, PubChemResolver
from opioidbench.chemistry.identifiers import normalize_cas, structure_formula_relation
from opioidbench.chemistry.structures import calculate_structure_features, scaffold_id


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _float_or_none(value: str | None) -> float | None:
    if value is None or not str(value).strip():
        return None
    try:
        return float(str(value).strip())
    except ValueError:
        return None


def _is_true(value: str | None) -> bool:
    return (value or "").strip().casefold() in {"1", "true", "yes"}


def _validate_pubchem_record(
    row: dict[str, str], record: EnrichmentRecord
) -> tuple[bool, str | None, str | None]:
    existing_inchikey = (row.get("inchikey") or "").strip().upper() or None
    source_conflict = any(
        _is_true(row.get(field))
        for field in (
            "name_conflict",
            "formula_conflict",
            "exact_mass_conflict",
            "structure_conflict",
        )
    )
    if source_conflict and not existing_inchikey:
        return False, None, "source_identity_conflict"

    source_formula = (row.get("molecular_formula") or "").strip() or None
    relation = None
    if source_formula and record.molecular_formula:
        relation = structure_formula_relation(source_formula, record.molecular_formula)
        if relation is None:
            return False, None, "source_formula_mismatch"

    if existing_inchikey and existing_inchikey != record.inchikey:
        return False, relation, "existing_structure_conflict"

    source_exact_mass = _float_or_none(row.get("source_monoisotopic_mass"))
    if source_exact_mass is not None and record.exact_mass is not None:
        if abs(source_exact_mass - record.exact_mass) > 0.0001 and relation == "exact_formula":
            return False, relation, "source_exact_mass_mismatch"
    return True, relation, None


def _empty_chebi_match(compound_id: str) -> ChEBIMatch:
    return ChEBIMatch(
        compound_id=compound_id,
        status="not_supplied",
        chebi_id=None,
        match_method=None,
        candidate_count=0,
        canonical_smiles=None,
        inchikey=None,
        formula_relation=None,
        error=None,
    )


def resolve_compound_sources(
    input_path: Path,
    *,
    output_dir: Path,
    chebi_sdf: Path | None = None,
    pubchem_cache_dir: Path | None = None,
    allow_pubchem_network: bool = False,
    max_records: int | None = None,
) -> dict[str, Any]:
    """Resolve ChEBI first, then PubChem CAS fallback, without guessing identities."""
    with input_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if max_records is not None:
        rows = rows[:max_records]

    output_dir.mkdir(parents=True, exist_ok=True)
    chebi_matches: list[ChEBIMatch]
    chebi_counts: dict[str, int] = {}
    if chebi_sdf is not None:
        if not chebi_sdf.is_file():
            raise FileNotFoundError(chebi_sdf)
        subset_input = output_dir / ".resolution_input.csv"
        fieldnames = list(rows[0].keys()) if rows else []
        with subset_input.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        try:
            chebi_matches, chebi_counts = resolve_chebi_compound_table(
                subset_input,
                chebi_sdf=chebi_sdf,
                output_path=output_dir / "chebi_matches.csv",
            )
        finally:
            subset_input.unlink(missing_ok=True)
    else:
        chebi_matches = [_empty_chebi_match(row.get("compound_id", "")) for row in rows]

    cache_root = pubchem_cache_dir or output_dir / "pubchem_cache"
    pubchem = PubChemResolver(
        MetadataCache(cache_root),
        allow_network=allow_pubchem_network,
    )

    resolved_rows: list[dict[str, Any]] = []
    review_rows: list[dict[str, Any]] = []
    descriptor_rows: list[dict[str, Any]] = []
    counts: dict[str, Any] = {
        "input_compounds": len(rows),
        "input_compounds_sha256": _sha256_file(input_path),
        "chebi_source_file": chebi_sdf.name if chebi_sdf is not None else None,
        "chebi_source_sha256": _sha256_file(chebi_sdf) if chebi_sdf is not None else None,
        "pubchem_network_enabled": allow_pubchem_network,
        "existing_structures": sum(bool(row.get("canonical_smiles")) for row in rows),
        "chebi_source_supplied": int(chebi_sdf is not None),
        "chebi_resolved": 0,
        "chebi_structure_filled": 0,
        "chebi_identity_only": 0,
        "chebi_ambiguous": 0,
        "chebi_conflict": 0,
        "pubchem_resolved": 0,
        "pubchem_structure_filled": 0,
        "pubchem_cache_miss_offline": 0,
        "pubchem_not_found": 0,
        "external_identity_unresolved": 0,
        "errors": 0,
    }

    for row, chebi_match in zip(rows, chebi_matches, strict=True):
        merged: dict[str, Any] = dict(row)
        existing_smiles = (row.get("canonical_smiles") or "").strip() or None
        existing_inchikey = (row.get("inchikey") or "").strip().upper() or None
        final_smiles = existing_smiles
        final_inchikey = existing_inchikey
        final_structure_source = (row.get("structure_source") or "").strip()
        final_formula_relation = (row.get("structure_formula_relation") or "").strip()
        chebi_id = ""
        pubchem_cid = (row.get("pubchem_cid") or "").strip()
        status = "unresolved"
        source = ""
        method = ""
        formula_relation = ""
        error = ""

        if chebi_match.status == "resolved":
            counts["chebi_resolved"] += 1
            chebi_id = chebi_match.chebi_id or ""
            source = "ChEBI"
            method = chebi_match.match_method or ""
            formula_relation = chebi_match.formula_relation or ""
            if existing_smiles:
                status = "chebi_identity_resolved_existing_structure"
            elif chebi_match.canonical_smiles:
                status = "chebi_structure_filled"
                counts["chebi_structure_filled"] += 1
                final_smiles = chebi_match.canonical_smiles
                final_inchikey = chebi_match.inchikey
                final_structure_source = "ChEBI bulk SDF"
                final_formula_relation = chebi_match.formula_relation or ""
            else:
                status = "chebi_identity_only"
                counts["chebi_identity_only"] += 1
        else:
            if chebi_match.status == "ambiguous":
                counts["chebi_ambiguous"] += 1
                error = chebi_match.error or "ambiguous_chebi_match"
            elif chebi_match.status == "conflict":
                counts["chebi_conflict"] += 1
                error = chebi_match.error or "chebi_consistency_conflict"

            cas = normalize_cas(row.get("cas"))
            if cas:
                try:
                    pubchem_record = pubchem.resolve_cas(cas)
                except (RuntimeError, ValueError) as exc:
                    pubchem_record = None
                    counts["errors"] += 1
                    error = f"{error} | PubChem:{exc}".strip(" |")
                if pubchem_record is not None:
                    ok, relation, validation_error = _validate_pubchem_record(row, pubchem_record)
                    if ok:
                        counts["pubchem_resolved"] += 1
                        pubchem_cid = pubchem_record.pubchem_cid
                        source = "PubChem PUG REST"
                        method = "cas"
                        formula_relation = relation or ""
                        if existing_smiles:
                            status = "pubchem_identity_resolved_existing_structure"
                        else:
                            status = "pubchem_structure_filled"
                            counts["pubchem_structure_filled"] += 1
                            final_smiles = pubchem_record.canonical_smiles
                            final_inchikey = pubchem_record.inchikey
                            final_structure_source = "PubChem PUG REST"
                            final_formula_relation = relation or ""
                    else:
                        error = f"{error} | PubChem:{validation_error}".strip(" |")
                elif allow_pubchem_network:
                    counts["pubchem_not_found"] += 1
                else:
                    counts["pubchem_cache_miss_offline"] += 1

        if status == "unresolved":
            counts["external_identity_unresolved"] += 1

        if final_smiles:
            features = calculate_structure_features(final_smiles)
            final_smiles = features.canonical_smiles
            if final_inchikey and final_inchikey != features.inchikey:
                counts["errors"] += 1
                error = f"{error} | final_structure_inchikey_conflict".strip(" |")
            final_inchikey = features.inchikey

            # External identity resolution must not silently rewrite scaffold labels
            # already accepted from a source-linked structure. Re-parsing an isomeric
            # SMILES can change RDKit's scaffold stereochemistry encoding even when the
            # molecular structure itself is unchanged. Preserve source-derived scaffold
            # fields when present, and calculate them only for newly filled structures
            # or legacy rows that lack scaffold metadata.
            existing_scaffold_smiles = (row.get("murcko_scaffold_smiles") or "").strip()
            existing_scaffold_id = (row.get("scaffold_id") or "").strip()
            if existing_smiles and existing_scaffold_id:
                final_scaffold_smiles = existing_scaffold_smiles
                final_scaffold = existing_scaffold_id
            else:
                final_scaffold_smiles = features.murcko_scaffold_smiles
                final_scaffold = scaffold_id(final_scaffold_smiles)

            descriptor_rows.append(
                {
                    "compound_id": row.get("compound_id", ""),
                    "canonical_smiles": features.canonical_smiles,
                    "inchikey": features.inchikey,
                    "molecular_weight_rdkit": features.molecular_weight,
                    "clogp": features.clogp,
                    "tpsa": features.tpsa,
                    "hbond_donors": features.hbond_donors,
                    "hbond_acceptors": features.hbond_acceptors,
                    "rotatable_bonds": features.rotatable_bonds,
                    "ring_count": features.ring_count,
                    "murcko_scaffold_smiles": final_scaffold_smiles,
                    "scaffold_id": final_scaffold,
                }
            )
            merged["canonical_smiles"] = final_smiles
            merged["inchikey"] = final_inchikey
            merged["murcko_scaffold_smiles"] = final_scaffold_smiles
            merged["scaffold_id"] = final_scaffold
            merged["structure_source"] = final_structure_source
            merged["structure_formula_relation"] = final_formula_relation

        merged.update(
            {
                "chebi_id": chebi_id,
                "pubchem_cid": pubchem_cid,
                "external_resolution_status": status,
                "external_resolution_source": source,
                "external_match_method": method,
                "external_formula_relation": formula_relation,
                "external_resolution_error": error,
            }
        )
        resolved_rows.append(merged)

        if status == "unresolved" or error or chebi_match.status in {"ambiguous", "conflict"}:
            review_rows.append(
                {
                    "compound_id": row.get("compound_id", ""),
                    "preferred_name": row.get("preferred_name", ""),
                    "cas": row.get("cas", ""),
                    "molecular_formula": row.get("molecular_formula", ""),
                    "inchikey": row.get("inchikey", ""),
                    "existing_structure": "yes" if existing_smiles else "no",
                    "chebi_status": chebi_match.status,
                    "chebi_candidate_count": chebi_match.candidate_count,
                    "resolution_status": status,
                    "error": error,
                    "recommended_action": (
                        "Supply the current ChEBI FULL SDF and/or a populated PubChem cache; "
                        "manually review ambiguity or conflicts rather than guessing."
                    ),
                }
            )

    existing_scaffold_ids = {
        (row.get("scaffold_id") or "").strip()
        for row in rows
        if (row.get("canonical_smiles") or "").strip() and (row.get("scaffold_id") or "").strip()
    }
    final_structured_rows = [
        row for row in resolved_rows if (row.get("canonical_smiles") or "").strip()
    ]
    final_scaffold_ids = {
        (row.get("scaffold_id") or "").strip()
        for row in final_structured_rows
        if (row.get("scaffold_id") or "").strip()
    }
    counts["final_structures"] = len(final_structured_rows)
    counts["final_structureless"] = len(rows) - len(final_structured_rows)
    counts["final_structure_coverage_fraction"] = (
        len(final_structured_rows) / len(rows) if rows else 0.0
    )
    counts["final_unique_scaffold_groups"] = len(final_scaffold_ids)
    counts["final_acyclic_structures"] = sum(
        (row.get("scaffold_id") or "").strip() == "murcko:acyclic" for row in final_structured_rows
    )
    counts["new_scaffold_groups"] = len(final_scaffold_ids - existing_scaffold_ids)

    conflict_reasons: Counter[str] = Counter()
    for match in chebi_matches:
        if not match.error:
            continue
        for item in match.error.split(" | "):
            conflict_reasons[item.rsplit(":", 1)[-1]] += 1
    counts["chebi_rejection_reason_counts"] = dict(sorted(conflict_reasons.items()))

    output_fields = list(rows[0].keys()) if rows else []
    for field in (
        "chebi_id",
        "pubchem_cid",
        "external_resolution_status",
        "external_resolution_source",
        "external_match_method",
        "external_formula_relation",
        "external_resolution_error",
    ):
        if field not in output_fields:
            output_fields.append(field)
    with (output_dir / "resolved_compounds.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(resolved_rows)

    descriptor_fields = [
        "compound_id",
        "canonical_smiles",
        "inchikey",
        "molecular_weight_rdkit",
        "clogp",
        "tpsa",
        "hbond_donors",
        "hbond_acceptors",
        "rotatable_bonds",
        "ring_count",
        "murcko_scaffold_smiles",
        "scaffold_id",
    ]
    with (output_dir / "calculated_descriptors.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=descriptor_fields)
        writer.writeheader()
        writer.writerows(descriptor_rows)

    review_fields = [
        "compound_id",
        "preferred_name",
        "cas",
        "molecular_formula",
        "inchikey",
        "existing_structure",
        "chebi_status",
        "chebi_candidate_count",
        "resolution_status",
        "error",
        "recommended_action",
    ]
    with (output_dir / "chemical_resolution_review.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=review_fields)
        writer.writeheader()
        writer.writerows(review_rows)

    if chebi_counts:
        counts["chebi_entries_scanned"] = chebi_counts["chebi_entries_scanned"]
    (output_dir / "resolution_report.json").write_text(
        json.dumps(counts, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return counts
