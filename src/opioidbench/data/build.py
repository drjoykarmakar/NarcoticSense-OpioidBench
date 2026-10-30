"""Versioned dataset construction and audit outputs."""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from dataclasses import replace
from datetime import date
from pathlib import Path
from typing import Any

from opioidbench import __version__
from opioidbench.chemistry.taxonomy import TaxonomyEntry, load_taxonomy
from opioidbench.data.download import sha256_file
from opioidbench.data.schema import NormalizedSpectrumRecord, RawSpectrum
from opioidbench.data.validation import QCReport, normalize_spectrum, validate_dataset
from opioidbench.experiments.reproducibility import write_json


def _taxonomy_label(entry: TaxonomyEntry | None) -> tuple[str | None, str | None]:
    if entry is None:
        return None, None
    status = entry.opioid_status.strip().lower()
    if entry.curation_status == "uncertain" or status == "uncertain":
        return "uncertain", entry.opioid_family or None
    if entry.curation_status != "curated":
        return None, None
    return entry.opioid_status or None, entry.opioid_family or None


def _apply_taxonomy(
    record: NormalizedSpectrumRecord,
    taxonomy: dict[str, TaxonomyEntry],
) -> NormalizedSpectrumRecord:
    label, family = _taxonomy_label(taxonomy.get(record.compound_id))
    return replace(record, opioid_label=label, opioid_family=family)


def _write_normalized(path: Path, records: list[NormalizedSpectrumRecord]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record.to_dict(), sort_keys=True) + "\n")


def _compound_rows(records: list[NormalizedSpectrumRecord]) -> list[dict[str, Any]]:
    groups: dict[str, list[NormalizedSpectrumRecord]] = defaultdict(list)
    for record in records:
        groups[record.compound_id].append(record)

    rows: list[dict[str, Any]] = []
    for compound_id in sorted(groups):
        group = groups[compound_id]
        raw_names = sorted(
            {
                name
                for item in group
                for name in (item.compound_name_raw, item.compound_name_companion)
                if name
            }
        )
        canonical_names = sorted(
            {item.compound_name_canonical for item in group if item.compound_name_canonical}
        )
        cas_values = sorted({item.cas for item in group if item.cas})
        formulas = sorted({item.molecular_formula for item in group if item.molecular_formula})
        exact_masses = sorted(
            {
                item.source_monoisotopic_mass
                for item in group
                if item.source_monoisotopic_mass is not None
            }
        )
        inchikeys = sorted({item.inchikey for item in group if item.inchikey})
        canonical_smiles = sorted(
            {item.canonical_smiles for item in group if item.canonical_smiles}
        )
        scaffold_ids = sorted({item.scaffold_id for item in group if item.scaffold_id})
        scaffold_smiles = sorted(
            {
                item.murcko_scaffold_smiles
                for item in group
                if item.murcko_scaffold_smiles is not None
            }
        )
        structure_sources = sorted(
            {item.structure_source for item in group if item.structure_source}
        )
        structure_relations = sorted(
            {item.structure_formula_relation for item in group if item.structure_formula_relation}
        )
        rows.append(
            {
                "compound_id": compound_id,
                "preferred_name": canonical_names[0] if len(canonical_names) == 1 else "",
                "aliases": " | ".join(raw_names),
                "cas": cas_values[0] if len(cas_values) == 1 else "",
                "molecular_formula": formulas[0] if len(formulas) == 1 else "",
                "source_monoisotopic_mass": (exact_masses[0] if len(exact_masses) == 1 else ""),
                "spectrum_count": len(group),
                "name_conflict": len(canonical_names) > 1,
                "formula_conflict": len(formulas) > 1,
                "exact_mass_conflict": len(exact_masses) > 1,
                "inchikey": inchikeys[0] if len(inchikeys) == 1 else "",
                "canonical_smiles": (canonical_smiles[0] if len(canonical_smiles) == 1 else ""),
                "scaffold_id": scaffold_ids[0] if len(scaffold_ids) == 1 else "",
                "murcko_scaffold_smiles": (scaffold_smiles[0] if len(scaffold_smiles) == 1 else ""),
                "structure_source": (structure_sources[0] if len(structure_sources) == 1 else ""),
                "structure_formula_relation": " | ".join(structure_relations),
                "structure_conflict": len(canonical_smiles) > 1,
                "structured_spectrum_count": sum(
                    item.canonical_smiles is not None for item in group
                ),
                "pubchem_cid": "",
            }
        )
    return rows


def _write_companion_review(path: Path, records: list[NormalizedSpectrumRecord]) -> None:
    """Write non-exact companion-name relationships for human review."""
    fieldnames = [
        "spectrum_id",
        "source_record_id",
        "compound_name_hpj",
        "compound_name_companion",
        "name_relation",
        "cas",
        "new_in_source_version",
    ]
    rows = [
        {
            "spectrum_id": record.spectrum_id,
            "source_record_id": record.source_record_id or "",
            "compound_name_hpj": record.compound_name_raw or "",
            "compound_name_companion": record.compound_name_companion or "",
            "name_relation": record.compound_name_resolution or "",
            "cas": record.cas or "",
            "new_in_source_version": record.new_in_source_version or "",
        }
        for record in records
        if record.compound_name_resolution and record.compound_name_resolution != "exact"
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_nist_structure_review(path: Path, rows: list[dict[str, Any]]) -> None:
    """Write local-only NIST structure mapping exceptions for human review."""
    fieldnames = [
        "source_record_id",
        "nist_library_id",
        "reason",
        "match_method",
        "source_formula",
        "structure_formula",
        "formula_relation",
        "source_exact_mass",
        "structure_exact_mass",
        "exact_mass_delta",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _write_compounds(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "compound_id",
        "preferred_name",
        "aliases",
        "cas",
        "molecular_formula",
        "source_monoisotopic_mass",
        "spectrum_count",
        "name_conflict",
        "formula_conflict",
        "exact_mass_conflict",
        "inchikey",
        "canonical_smiles",
        "scaffold_id",
        "murcko_scaffold_smiles",
        "structure_source",
        "structure_formula_relation",
        "structure_conflict",
        "structured_spectrum_count",
        "pubchem_cid",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _taxonomy_counts(taxonomy: dict[str, TaxonomyEntry]) -> dict[str, Any]:
    curated = [entry for entry in taxonomy.values() if entry.curation_status == "curated"]
    uncertain = [entry for entry in taxonomy.values() if entry.curation_status == "uncertain"]
    opioid = [entry for entry in curated if entry.opioid_status.strip().lower() == "opioid"]
    nonopioid = [
        entry
        for entry in curated
        if entry.opioid_status.strip().lower() in {"non-opioid", "nonopioid"}
    ]
    families = {entry.opioid_family for entry in opioid if entry.opioid_family}
    return {
        "number_taxonomy_entries": len(taxonomy),
        "number_curated_taxonomy_entries": len(curated),
        "number_uncertain_labels": len(uncertain),
        "number_opioid_labeled_compounds": len(opioid),
        "number_nonopioid_labeled_compounds": len(nonopioid),
        "number_opioid_families": len(families),
    }


def _audit_payload(
    records: list[NormalizedSpectrumRecord],
    report: QCReport,
    taxonomy: dict[str, TaxonomyEntry],
) -> dict[str, Any]:
    compound_ids = {record.compound_id for record in records}
    compounds_with_structures = {
        record.compound_id for record in records if record.canonical_smiles and record.inchikey
    }
    compound_structure_sets: dict[str, set[str]] = defaultdict(set)
    compound_scaffold_sets: dict[str, set[str]] = defaultdict(set)
    for record in records:
        if record.canonical_smiles and record.inchikey:
            compound_structure_sets[record.compound_id].add(record.canonical_smiles)
            if record.scaffold_id:
                compound_scaffold_sets[record.compound_id].add(record.scaffold_id)
    unambiguous_structure_compounds = {
        compound_id
        for compound_id, structures in compound_structure_sets.items()
        if len(structures) == 1
    }
    conflicting_structure_compounds = {
        compound_id
        for compound_id, structures in compound_structure_sets.items()
        if len(structures) > 1
    }
    unambiguous_compound_scaffolds = {
        scaffold
        for compound_id in unambiguous_structure_compounds
        for scaffold in compound_scaffold_sets.get(compound_id, set())
    }
    spectra_with_structures = [record for record in records if record.canonical_smiles]
    structure_sources = Counter(
        record.structure_source for record in spectra_with_structures if record.structure_source
    )
    structure_methods = Counter(
        record.structure_match_method
        for record in spectra_with_structures
        if record.structure_match_method
    )
    structure_relations = Counter(
        record.structure_formula_relation
        for record in spectra_with_structures
        if record.structure_formula_relation
    )
    scaffold_ids = {record.scaffold_id for record in spectra_with_structures if record.scaffold_id}
    taxonomy_counts = _taxonomy_counts(taxonomy)
    companion_relations = Counter(
        record.compound_name_resolution for record in records if record.compound_name_resolution
    )
    new_source_versions = Counter(
        record.new_in_source_version for record in records if record.new_in_source_version
    )
    exact_mass_records = [
        record for record in records if record.source_monoisotopic_mass is not None
    ]
    source_comment_relations = Counter(
        record.instrument_metadata.get("SOURCE_COMMENT_RELATION")
        for record in records
        if record.instrument_metadata.get("SOURCE_COMMENT_RELATION")
    )
    source_comment_full_records = sum(
        bool(record.instrument_metadata.get("SOURCE_COMMENT_FULL")) for record in records
    )
    nominal_mass_mismatches = [
        record
        for record in exact_mass_records
        if record.molecular_weight is not None
        and round(record.source_monoisotopic_mass) != record.molecular_weight
    ]
    return {
        "number_spectra": len(records),
        "number_unique_compounds": len(compound_ids),
        "number_compounds_with_structures": len(compounds_with_structures),
        "number_compounds_with_unambiguous_structures": len(unambiguous_structure_compounds),
        "number_compounds_with_structure_conflicts": len(conflicting_structure_compounds),
        "number_spectra_with_structures": len(spectra_with_structures),
        "number_unique_murcko_scaffolds": len(scaffold_ids),
        "number_unique_unambiguous_compound_scaffolds": len(unambiguous_compound_scaffolds),
        "number_acyclic_structured_spectra": sum(
            record.scaffold_id == "murcko:acyclic" for record in spectra_with_structures
        ),
        "structure_source_counts": dict(sorted(structure_sources.items())),
        "structure_match_method_counts": dict(sorted(structure_methods.items())),
        "structure_formula_relation_counts": dict(sorted(structure_relations.items())),
        "number_invalid_spectra": report.number_invalid,
        "number_valid_spectra": report.number_valid,
        "placeholder_cas_records": report.issue_counts.get("placeholder_cas", 0),
        "invalid_cas_records": report.issue_counts.get("invalid_cas", 0),
        "duplicate_spectrum_groups": report.issue_counts.get("duplicate_spectrum", 0),
        "duplicate_chemical_identity_groups": report.issue_counts.get(
            "duplicate_chemical_identity", 0
        ),
        "conflicting_compound_name_groups": report.issue_counts.get(
            "conflicting_compound_names", 0
        ),
        "companion_metadata_records": sum(companion_relations.values()),
        "companion_name_relation_counts": dict(sorted(companion_relations.items())),
        "new_in_source_version_counts": dict(sorted(new_source_versions.items())),
        "source_monoisotopic_mass_records": len(exact_mass_records),
        "source_m_plus_proton_mass_records": sum(
            record.source_m_plus_proton_mass is not None for record in records
        ),
        "source_comment_full_records": source_comment_full_records,
        "source_comment_relation_counts": dict(sorted(source_comment_relations.items())),
        "source_comment_extended_records": source_comment_relations.get("hpj_truncated", 0),
        "exact_mass_nominal_weight_rounding_mismatch_records": len(nominal_mass_mismatches),
        **taxonomy_counts,
        "taxonomy_complete": False,
    }


def _write_qc_markdown(path: Path, report: QCReport, audit: dict[str, Any]) -> None:
    lines = [
        "# Dataset QC summary",
        "",
        f"- Spectra: {audit['number_spectra']}",
        f"- Valid spectra: {audit['number_valid_spectra']}",
        f"- Invalid spectra: {audit['number_invalid_spectra']}",
        f"- Conservative unique compounds: {audit['number_unique_compounds']}",
        f"- Compounds with structure evidence: {audit['number_compounds_with_structures']}",
        (
            "- Compounds with one unambiguous structure: "
            f"{audit['number_compounds_with_unambiguous_structures']}"
        ),
        (
            "- Compound identities with conflicting structures: "
            f"{audit['number_compounds_with_structure_conflicts']}"
        ),
        f"- Spectra with structures: {audit['number_spectra_with_structures']}",
        f"- Spectrum-level unique Murcko scaffold IDs: {audit['number_unique_murcko_scaffolds']}",
        (
            "- Unambiguous compound-level Murcko scaffold IDs: "
            f"{audit['number_unique_unambiguous_compound_scaffolds']}"
        ),
        f"- Structure formula relations: {audit['structure_formula_relation_counts']}",
        f"- Placeholder CAS records: {audit['placeholder_cas_records']}",
        f"- Invalid CAS records: {audit['invalid_cas_records']}",
        f"- Exact duplicate spectrum groups: {audit['duplicate_spectrum_groups']}",
        f"- Repeated chemical-identity groups: {audit['duplicate_chemical_identity_groups']}",
        f"- Conflicting-name identity groups: {audit['conflicting_compound_name_groups']}",
        f"- Companion metadata records: {audit['companion_metadata_records']}",
        f"- Companion name relations: {audit['companion_name_relation_counts']}",
        f"- New-in-version records: {audit['new_in_source_version_counts']}",
        f"- Source monoisotopic-mass records: {audit['source_monoisotopic_mass_records']}",
        f"- Source M+proton records: {audit['source_m_plus_proton_mass_records']}",
        f"- Full source-comment records: {audit['source_comment_full_records']}",
        f"- Source-comment relations: {audit['source_comment_relation_counts']}",
        f"- HPJ source-comments extended from Agilent: {audit['source_comment_extended_records']}",
        (
            "- Exact-mass vs nominal-MW rounding mismatches: "
            f"{audit['exact_mass_nominal_weight_rounding_mismatch_records']}"
        ),
        f"- Curated taxonomy entries: {audit['number_curated_taxonomy_entries']}",
        f"- Uncertain taxonomy entries: {audit['number_uncertain_labels']}",
        "",
        "## QC issue counts",
        "",
    ]
    if report.issue_counts:
        lines.extend(f"- `{code}`: {count}" for code, count in report.issue_counts.items())
    else:
        lines.append("- No QC issues detected.")
    lines.extend(
        [
            "",
            (
                "> Taxonomy coverage is incomplete. Unlabeled compounds must not be "
                "treated as non-opioids."
            ),
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_dataset(
    spectra: list[RawSpectrum],
    *,
    output_dir: Path,
    dataset_version: str,
    source_file: Path,
    source_file_sha256: str,
    taxonomy_path: Path | None = None,
    source_archive: Path | None = None,
    source_archive_sha256: str | None = None,
    companion_file: Path | None = None,
    companion_file_sha256: str | None = None,
    companion_alignment: dict[str, Any] | None = None,
    exact_mass_file: Path | None = None,
    exact_mass_file_sha256: str | None = None,
    exact_mass_alignment: dict[str, Any] | None = None,
    agilent_archive: Path | None = None,
    agilent_archive_sha256: str | None = None,
    agilent_header_alignment: dict[str, Any] | None = None,
    nist_archive: Path | None = None,
    nist_archive_sha256: str | None = None,
    nist_structure_member: str | None = None,
    nist_structure_member_sha256: str | None = None,
    nist_structure_alignment: dict[str, Any] | None = None,
    nist_structure_review: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Validate, normalize, inventory, audit, and version one dataset build."""
    taxonomy = load_taxonomy(taxonomy_path) if taxonomy_path else {}
    validated, report = validate_dataset(spectra)
    normalized = [_apply_taxonomy(normalize_spectrum(item), taxonomy) for item in validated]

    output_dir.mkdir(parents=True, exist_ok=True)
    normalized_path = output_dir / "normalized_records.jsonl"
    qc_path = output_dir / "data_qc.json"
    audit_path = output_dir / "dataset_audit.json"
    compounds_path = output_dir / "compounds.csv"
    manifest_path = output_dir / "dataset_manifest.json"
    qc_md_path = output_dir / "data_qc.md"
    companion_review_path = output_dir / "companion_name_review.csv"
    nist_review_path = output_dir / "nist_structure_review.csv"

    _write_normalized(normalized_path, normalized)
    write_json(qc_path, {**report.to_dict(), "source_file_sha256": source_file_sha256})
    compound_rows = _compound_rows(normalized)
    _write_compounds(compounds_path, compound_rows)
    if companion_file:
        _write_companion_review(companion_review_path, normalized)
    if nist_archive is not None:
        _write_nist_structure_review(nist_review_path, nist_structure_review or [])
    audit = _audit_payload(normalized, report, taxonomy)
    if nist_structure_alignment is not None:
        audit["nist_structure_alignment"] = nist_structure_alignment
    write_json(audit_path, audit)
    _write_qc_markdown(qc_md_path, report, audit)

    manifest = {
        "dataset_version": dataset_version,
        "creation_date": date.today().isoformat(),
        "source_versions": sorted(
            {f"{record.source}:{record.source_version or 'unknown'}" for record in normalized}
        ),
        "source_file": source_file.name,
        "source_file_sha256": source_file_sha256,
        "source_archive": source_archive.name if source_archive else None,
        "source_archive_sha256": source_archive_sha256,
        "companion_file": companion_file.name if companion_file else None,
        "companion_file_sha256": companion_file_sha256,
        "companion_alignment": companion_alignment,
        "exact_mass_file": exact_mass_file.name if exact_mass_file else None,
        "exact_mass_file_sha256": exact_mass_file_sha256,
        "exact_mass_alignment": exact_mass_alignment,
        "agilent_archive": agilent_archive.name if agilent_archive else None,
        "agilent_archive_sha256": agilent_archive_sha256,
        "agilent_header_alignment": agilent_header_alignment,
        "nist_archive": nist_archive.name if nist_archive else None,
        "nist_archive_sha256": nist_archive_sha256,
        "nist_structure_member": nist_structure_member,
        "nist_structure_member_sha256": nist_structure_member_sha256,
        "nist_structure_alignment": nist_structure_alignment,
        "nist_structure_interpretation": (
            "same SWGDRUG release in NIST format; linked molfiles are accepted only "
            "after conservative ID/name/CAS plus formula consistency checks"
            if nist_archive
            else None
        ),
        "agilent_header_interpretation": (
            "same SWGDRUG release in vendor format; used only to restore full source comments"
            if agilent_archive
            else None
        ),
        "exact_mass_interpretation": (
            "source-provided formula-derived Agilent PCDL value; workbook states unverified"
            if exact_mass_file
            else None
        ),
        "number_spectra": audit["number_spectra"],
        "number_compounds": audit["number_unique_compounds"],
        "number_opioids": audit["number_opioid_labeled_compounds"],
        "number_nonopioids": audit["number_nonopioid_labeled_compounds"],
        "number_families": audit["number_opioid_families"],
        "taxonomy_complete": audit["taxonomy_complete"],
        "hash": sha256_file(normalized_path),
        "pipeline_version": __version__,
        "outputs": {
            "normalized_records": normalized_path.name,
            "quality_control": qc_path.name,
            "quality_control_summary": qc_md_path.name,
            "compound_inventory": compounds_path.name,
            "dataset_audit": audit_path.name,
            "companion_name_review": (companion_review_path.name if companion_file else None),
            "nist_structure_review": (nist_review_path.name if nist_archive else None),
        },
    }
    write_json(manifest_path, manifest)
    return {"audit": audit, "manifest": manifest, "report": report.to_dict()}
