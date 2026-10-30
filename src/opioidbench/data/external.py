"""Conservative construction of ENFSI external-test candidates.

ENFSI source data are redistribution-restricted. Public outputs from this module are
aggregate-only; the optional row-level manifest is explicitly local/private.
"""

from __future__ import annotations

import csv
import json
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from opioidbench.chemistry.chebi import collect_relevant_chebi_entries, select_chebi_match
from opioidbench.chemistry.evidence import OpioidEvidence, load_opioid_evidence
from opioidbench.chemistry.identifiers import (
    formula_atom_counts,
    normalize_cas,
    normalize_compound_name,
)
from opioidbench.chemistry.structures import calculate_structure_features, scaffold_id
from opioidbench.data.enfsi import read_enfsi_entry_list
from opioidbench.data.overlap import audit_enfsi_spectral_overlap
from opioidbench.data.schema import RawSpectrum


@dataclass(frozen=True, slots=True)
class SupplementalStructure:
    preferred_name: str
    aliases: tuple[str, ...]
    molecular_formula: str
    canonical_smiles: str
    inchikey: str
    public_id: str
    source: str
    source_url: str

    @property
    def name_keys(self) -> frozenset[str]:
        keys = {_name_key(value) for value in (self.preferred_name, *self.aliases)}
        return frozenset(key for key in keys if key)


def _name_key(value: str | None) -> str | None:
    normalized = normalize_compound_name(value)
    return normalized.casefold() if normalized else None


def _formula_key(value: str | None) -> tuple[tuple[str, int], ...] | None:
    counts = formula_atom_counts(value)
    if counts is None:
        return None
    return tuple(sorted(counts.items()))


def _split_pipe(value: str | None) -> tuple[str, ...]:
    if not value:
        return ()
    return tuple(part.strip() for part in value.split("|") if part.strip())


def _aliases(name: str) -> tuple[str, ...]:
    values: list[str] = [name]
    values.extend(part.strip() for part in name.split(";") if part.strip())
    values.extend(
        match.group(1).strip()
        for match in re.finditer(r"\(([^()]*)\)", name)
        if match.group(1).strip()
    )
    outside = re.sub(r"\([^()]*\)", " ", name).strip()
    if outside:
        values.append(outside)
    return tuple(dict.fromkeys(values))


def _enfsi_resolution_rows(
    spectra: list[RawSpectrum], entry_list_path: Path
) -> list[dict[str, str]]:
    entries = read_enfsi_entry_list(entry_list_path)
    if len(entries) != len(spectra):
        raise ValueError("ENFSI entry list and spectra differ in count")
    rows: list[dict[str, str]] = []
    for spectrum, entry in zip(spectra, entries, strict=True):
        expected_id = f"Library Entry {entry.source_index}"
        if spectrum.provenance.source_record_id != expected_id:
            raise ValueError("ENFSI entry list and spectra do not align by source record ID")
        aliases = _aliases(entry.name)
        rows.append(
            {
                "compound_id": f"enfsi:{entry.source_index}",
                "preferred_name": entry.name,
                "aliases": "|".join(aliases),
                "cas": normalize_cas(spectrum.metadata.get("CAS_REGISTRY_NO")) or "",
                "molecular_formula": entry.formula or "",
                "source_monoisotopic_mass": "",
                "inchikey": "",
                "canonical_smiles": "",
                "name_conflict": "False",
                "formula_conflict": "False",
                "exact_mass_conflict": "False",
                "structure_conflict": "False",
            }
        )
    return rows


def _load_swgdrug_structure_sets(path: Path) -> tuple[set[str], set[str], int]:
    inchikeys: set[str] = set()
    scaffolds: set[str] = set()
    usable_rows = 0
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if (row.get("structure_conflict") or "").strip().casefold() == "true":
                continue
            inchikey = (row.get("inchikey") or "").strip().upper()
            scaffold = (row.get("scaffold_id") or "").strip()
            if inchikey or scaffold:
                usable_rows += 1
            if inchikey:
                inchikeys.add(inchikey)
            if scaffold:
                scaffolds.add(scaffold)
    return inchikeys, scaffolds, usable_rows


def _peak_signature(spectrum: RawSpectrum) -> tuple[tuple[float, float], ...]:
    return tuple((float(peak.mz), float(peak.intensity)) for peak in spectrum.peaks)


def _exact_duplicate_indices(spectra: list[RawSpectrum]) -> set[int]:
    groups: dict[tuple[tuple[float, float], ...], list[int]] = defaultdict(list)
    for idx, spectrum in enumerate(spectra):
        groups[_peak_signature(spectrum)].append(idx)
    duplicates: set[int] = set()
    for indices in groups.values():
        if len(indices) > 1:
            duplicates.update(indices)
    return duplicates


def _structure_formula(smiles: str) -> str:
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise RuntimeError("RDKit is required for external structure validation") from exc
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid supplemental SMILES: {smiles!r}")
    return str(rdMolDescriptors.CalcMolFormula(mol))


def _load_supplemental_structures(
    path: Path | None,
) -> tuple[dict[tuple[str, tuple[tuple[str, int], ...]], list[SupplementalStructure]], int]:
    """Load public structures indexed only by exact normalized name/alias and formula."""
    if path is None:
        return {}, 0
    index: dict[tuple[str, tuple[tuple[str, int], ...]], list[SupplementalStructure]] = defaultdict(
        list
    )
    count = 0
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "preferred_name",
            "aliases",
            "molecular_formula",
            "canonical_smiles",
            "inchikey",
            "pubchem_cid",
            "source",
            "source_url",
        }
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValueError(
                "Supplemental structure file is missing columns: " + ", ".join(missing)
            )
        for raw in reader:
            preferred_name = (raw.get("preferred_name") or "").strip()
            formula = (raw.get("molecular_formula") or "").strip()
            smiles = (raw.get("canonical_smiles") or "").strip()
            supplied_inchikey = (raw.get("inchikey") or "").strip().upper()
            if not preferred_name or not formula or not smiles or not supplied_inchikey:
                raise ValueError(
                    "Supplemental structure rows require name, formula, SMILES, InChIKey"
                )
            features = calculate_structure_features(smiles)
            if features.inchikey.upper() != supplied_inchikey:
                raise ValueError(f"Supplemental structure InChIKey mismatch for {preferred_name!r}")
            source_formula = _formula_key(formula)
            calculated_formula = _formula_key(_structure_formula(features.canonical_smiles))
            if source_formula is None or source_formula != calculated_formula:
                raise ValueError(f"Supplemental structure formula mismatch for {preferred_name!r}")
            item = SupplementalStructure(
                preferred_name=preferred_name,
                aliases=_split_pipe(raw.get("aliases")),
                molecular_formula=formula,
                canonical_smiles=features.canonical_smiles,
                inchikey=features.inchikey.upper(),
                public_id=(raw.get("pubchem_cid") or "").strip(),
                source=(raw.get("source") or "").strip(),
                source_url=(raw.get("source_url") or "").strip(),
            )
            formula_key = _formula_key(formula)
            assert formula_key is not None
            for name_key in item.name_keys:
                index[(name_key, formula_key)].append(item)
            count += 1
    return dict(index), count


def _supplemental_structure_for_row(
    row: dict[str, str],
    index: dict[tuple[str, tuple[tuple[str, int], ...]], list[SupplementalStructure]],
) -> SupplementalStructure | None:
    formula_key = _formula_key(row.get("molecular_formula"))
    if formula_key is None:
        return None
    names = {
        _name_key(value)
        for value in (row.get("preferred_name", ""), *_split_pipe(row.get("aliases")))
    }
    matches: dict[str, SupplementalStructure] = {}
    for name_key in names:
        if not name_key:
            continue
        for item in index.get((name_key, formula_key), []):
            matches[item.inchikey] = item
    if len(matches) != 1:
        return None
    return next(iter(matches.values()))


def _label_row_from_evidence(
    row: dict[str, str], evidence: list[OpioidEvidence]
) -> tuple[str, tuple[str, ...], tuple[str, ...]]:
    """Return exact evidence label, matched evidence IDs, and match methods."""
    cas = normalize_cas(row.get("cas"))
    names = {
        _name_key(value)
        for value in (row.get("preferred_name", ""), *_split_pipe(row.get("aliases")))
    }
    matches: list[tuple[OpioidEvidence, str]] = []
    for item in evidence:
        if cas and item.cas and cas == item.cas:
            matches.append((item, "cas"))
            continue
        if any(name and name in item.name_keys for name in names):
            matches.append((item, "name"))
    if not matches:
        return "pending", (), ()
    statuses = {item.opioid_status for item, _ in matches}
    label = "conflict" if len(statuses) != 1 else next(iter(statuses))
    evidence_ids = tuple(sorted({item.evidence_id for item, _ in matches}))
    methods = tuple(sorted({method for _, method in matches}))
    return label, evidence_ids, methods


def audit_enfsi_external_candidates(
    enfsi_spectra: list[RawSpectrum],
    swgdrug_spectra: list[RawSpectrum],
    *,
    entry_list_path: Path,
    normalized_swgdrug_path: Path,
    resolved_swgdrug_compounds_path: Path,
    chebi_sdf_path: Path,
    supplemental_structures_path: Path | None = None,
    taxonomy_evidence_path: Path | None = None,
    near_identical_threshold: float = 0.9999,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Build an aggregate audit and a restricted local candidate manifest.

    Candidate structures come first from conservative ChEBI resolution. A public
    supplemental structure may fill a ChEBI-unresolved record only after an exact
    normalized name/alias plus molecular-formula match. CAS alone never selects a
    supplemental structure.
    """
    overlap_report, spectral_details = audit_enfsi_spectral_overlap(
        enfsi_spectra,
        swgdrug_spectra,
        near_identical_threshold=near_identical_threshold,
        entry_list_path=entry_list_path,
        normalized_swgdrug_path=normalized_swgdrug_path,
    )
    rows = _enfsi_resolution_rows(enfsi_spectra, entry_list_path)
    candidates, scanned = collect_relevant_chebi_entries(chebi_sdf_path, rows)
    matches = [select_chebi_match(row, candidates.get(idx, [])) for idx, row in enumerate(rows)]
    supplemental_index, supplemental_records = _load_supplemental_structures(
        supplemental_structures_path
    )
    evidence = load_opioid_evidence(taxonomy_evidence_path) if taxonomy_evidence_path else []
    swg_inchikeys, swg_scaffolds, swg_structured_rows = _load_swgdrug_structure_sets(
        resolved_swgdrug_compounds_path
    )
    duplicate_indices = _exact_duplicate_indices(enfsi_spectra)

    private_rows: list[dict[str, Any]] = []
    candidate_indices: list[int] = []
    status_counts: Counter[str] = Counter(match.status for match in matches)
    base_count = 0
    base_chebi_resolved = 0
    base_supplemental_resolved = 0
    supplemental_used_all = 0
    base_structure_resolved = 0
    base_identity_overlap = 0
    base_scaffold_overlap = 0
    base_internal_duplicate = 0

    for idx, (row, detail, match) in enumerate(zip(rows, spectral_details, matches, strict=True)):
        base = detail.metadata_overlap_status == "unmatched" and not detail.near_identical
        if base:
            base_count += 1
        canonical_smiles = match.canonical_smiles if match.status == "resolved" else None
        inchikey = (match.inchikey or "").upper() if match.status == "resolved" else ""
        structure_source = "ChEBI" if canonical_smiles else ""
        public_structure_id: str | None = None
        if not canonical_smiles and supplemental_index:
            supplemental = _supplemental_structure_for_row(row, supplemental_index)
            if supplemental is not None:
                canonical_smiles = supplemental.canonical_smiles
                inchikey = supplemental.inchikey
                structure_source = supplemental.source or "public_supplement"
                public_structure_id = supplemental.public_id or None
                supplemental_used_all += 1
                if base:
                    base_supplemental_resolved += 1

        scaffold = ""
        if canonical_smiles:
            features = calculate_structure_features(canonical_smiles)
            canonical_smiles = features.canonical_smiles
            inchikey = features.inchikey.upper()
            scaffold = scaffold_id(features.murcko_scaffold_smiles)

        if base and match.status == "resolved":
            base_chebi_resolved += 1
        if base and canonical_smiles and inchikey and scaffold:
            base_structure_resolved += 1

        identity_overlap = bool(inchikey and inchikey in swg_inchikeys)
        scaffold_overlap = bool(scaffold and scaffold in swg_scaffolds)
        internal_duplicate = idx in duplicate_indices
        if base and identity_overlap:
            base_identity_overlap += 1
        if base and scaffold_overlap:
            base_scaffold_overlap += 1
        if base and internal_duplicate:
            base_internal_duplicate += 1

        evidence_label, evidence_ids, evidence_methods = _label_row_from_evidence(row, evidence)
        candidate = bool(
            base
            and canonical_smiles
            and inchikey
            and scaffold
            and not identity_overlap
            and not scaffold_overlap
            and not internal_duplicate
        )
        if candidate:
            candidate_indices.append(idx)

        private_rows.append(
            {
                "enfsi_source_record_id": detail.enfsi_source_record_id,
                "best_swgdrug_cosine": detail.cosine_similarity,
                "metadata_overlap_status": detail.metadata_overlap_status,
                "near_identical": detail.near_identical,
                "chebi_status": match.status,
                "chebi_id": match.chebi_id,
                "chebi_match_method": match.match_method,
                "structure_resolution_source": structure_source or None,
                "public_structure_id": public_structure_id,
                "inchikey": inchikey or None,
                "scaffold_id": scaffold or None,
                "identity_overlap_swgdrug": identity_overlap,
                "scaffold_overlap_swgdrug": scaffold_overlap,
                "enfsi_exact_duplicate_group": internal_duplicate,
                "taxonomy_label": evidence_label,
                "taxonomy_evidence_ids": list(evidence_ids),
                "taxonomy_match_methods": list(evidence_methods),
                "external_candidate_record": candidate,
                "primary_unique_identity_representative": False,
            }
        )

    first_by_inchikey: dict[str, int] = {}
    for idx in candidate_indices:
        inchikey = str(private_rows[idx]["inchikey"])
        first_by_inchikey.setdefault(inchikey, idx)
    for idx in first_by_inchikey.values():
        private_rows[idx]["primary_unique_identity_representative"] = True

    candidate_scaffolds = {
        str(private_rows[idx]["scaffold_id"])
        for idx in candidate_indices
        if private_rows[idx]["scaffold_id"]
    }
    candidate_label_counts: Counter[str] = Counter(
        str(private_rows[idx]["taxonomy_label"]) for idx in first_by_inchikey.values()
    )
    report: dict[str, Any] = {
        "enfsi_entries": len(enfsi_spectra),
        "swgdrug_spectra": len(swgdrug_spectra),
        "near_identical_threshold": near_identical_threshold,
        "metadata_unmatched_spectrum_disjoint_entries": base_count,
        "chebi_entries_scanned": scanned,
        "chebi_status_counts_all_enfsi": dict(sorted(status_counts.items())),
        "supplemental_structure_records_loaded": supplemental_records,
        "supplemental_structures_used_all_enfsi": supplemental_used_all,
        "base_candidates_chebi_resolved": base_chebi_resolved,
        "base_candidates_supplemental_structure_resolved": base_supplemental_resolved,
        "base_candidates_structure_resolved": base_structure_resolved,
        "base_candidates_exact_inchikey_overlap_swgdrug": base_identity_overlap,
        "base_candidates_exact_scaffold_overlap_swgdrug": base_scaffold_overlap,
        "base_candidates_enfsi_exact_duplicate_group_records": base_internal_duplicate,
        "external_candidate_records": len(candidate_indices),
        "external_candidate_unique_inchikeys": len(first_by_inchikey),
        "external_candidate_unique_scaffolds": len(candidate_scaffolds),
        "primary_unique_identity_representatives": len(first_by_inchikey),
        "external_candidate_labeled_opioid": candidate_label_counts.get("opioid", 0),
        "external_candidate_labeled_nonopioid": candidate_label_counts.get("non-opioid", 0),
        "external_candidate_label_pending": candidate_label_counts.get("pending", 0),
        "external_candidate_label_conflict": candidate_label_counts.get("conflict", 0),
        "swgdrug_unambiguous_structure_rows": swg_structured_rows,
        "swgdrug_unique_inchikeys": len(swg_inchikeys),
        "swgdrug_unique_scaffolds": len(swg_scaffolds),
        "restricted_source_data_in_aggregate_output": False,
        "candidate_definition": (
            "ENFSI metadata-unmatched to SWGDRUG; best SWGDRUG cosine below threshold; "
            "structure resolved by conservative ChEBI matching or exact-name-plus-formula "
            "public supplement; exact InChIKey absent from SWGDRUG; exact Murcko scaffold "
            "absent from SWGDRUG; not in an ENFSI exact duplicate-spectrum group."
        ),
        "interpretation": (
            "These are conservative external-test candidates, not proof of independent "
            "acquisition or complete chemical novelty. Primary evaluation uses one "
            "representative per exact InChIKey. Only source-backed labels are counted; "
            "unmatched candidates remain pending rather than being treated as negatives."
        ),
        "spectral_overlap_summary": {
            "near_identical_best_match_records": overlap_report.near_identical_best_match_records,
            "unmatched_metadata_entries": overlap_report.unmatched_metadata_entries,
            "metadata_unmatched_spectrum_disjoint_candidates": (
                overlap_report.metadata_unmatched_spectrum_disjoint_candidates
            ),
        },
    }
    return report, private_rows


def write_private_external_candidate_manifest(path: Path, rows: list[dict[str, Any]]) -> None:
    """Write restricted local row-level candidate decisions; never a public artifact."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")
