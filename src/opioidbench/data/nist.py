"""Conservative extraction of SWGDRUG structures from the NIST-format archive."""

from __future__ import annotations

import hashlib
import re
import zipfile
from collections import Counter, defaultdict
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

from opioidbench.chemistry.identifiers import (
    normalize_cas,
    normalize_compound_name,
    structure_formula_relation,
)
from opioidbench.chemistry.structures import ChemistryDependencyError, scaffold_id
from opioidbench.data.schema import RawSpectrum

_LIBRARY_ID_RE = re.compile(r"Library ID\s*=\s*(\d+)", re.IGNORECASE)
_CAS_RN_RE = re.compile(r"CAS rn\s*=\s*(\d+)", re.IGNORECASE)


@dataclass(frozen=True, slots=True)
class NistStructureEntry:
    """One library-linked structure recovered from ``USRSTRUC.DB``."""

    library_id: int
    source_name: str
    source_cas: str | None
    canonical_smiles: str | None
    inchikey: str | None
    molecular_formula: str | None
    exact_mass: float | None
    murcko_scaffold_smiles: str | None
    scaffold_id: str | None
    parse_error: str | None


@dataclass(frozen=True, slots=True)
class NistStructureAlignment:
    """Aggregate, redistribution-safe structure alignment statistics."""

    linked_structure_blocks: int
    parsed_structure_blocks: int
    parse_failures: int
    matched_spectra: int
    unmatched_spectra: int
    conflicting_spectra: int
    duplicate_candidate_collapses: int
    match_method_counts: dict[str, int]
    formula_relation_counts: dict[str, int]
    exact_mass_checked_records: int
    exact_mass_warning_records: int

    def to_dict(self) -> dict[str, Any]:
        return {
            "linked_structure_blocks": self.linked_structure_blocks,
            "parsed_structure_blocks": self.parsed_structure_blocks,
            "parse_failures": self.parse_failures,
            "matched_spectra": self.matched_spectra,
            "unmatched_spectra": self.unmatched_spectra,
            "conflicting_spectra": self.conflicting_spectra,
            "duplicate_candidate_collapses": self.duplicate_candidate_collapses,
            "match_method_counts": self.match_method_counts,
            "formula_relation_counts": self.formula_relation_counts,
            "exact_mass_checked_records": self.exact_mass_checked_records,
            "exact_mass_warning_records": self.exact_mass_warning_records,
        }


def _rdkit():
    try:
        from rdkit import Chem, rdBase
        from rdkit.Chem import Descriptors, rdMolDescriptors
        from rdkit.Chem.Scaffolds import MurckoScaffold
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise ChemistryDependencyError(
            "RDKit is required for NIST structure extraction; install opioidbench[chem]."
        ) from exc
    return Chem, rdBase, Descriptors, rdMolDescriptors, MurckoScaffold


def _cas_from_digits(value: str) -> str | None:
    digits = re.sub(r"\D", "", value)
    if len(digits) < 5:
        return None
    candidate = f"{int(digits[:-3])}-{digits[-3:-1]}-{digits[-1]}"
    return normalize_cas(candidate)


def _molblock(block: str) -> str:
    """Convert the older NIST molfile variant to V2000 syntax when needed."""
    if "V2000" in block:
        return block + "\n"
    lines = block.splitlines()
    if len(lines) < 4:
        return block + "\n"
    lines[3] = lines[3].rstrip() + " V2000"
    if not any(line.strip() == "M  END" for line in lines):
        lines.append("M  END")
    return "\n".join(lines) + "\n"


def _parse_structure_block(block: str) -> NistStructureEntry | None:
    match = _LIBRARY_ID_RE.search(block)
    if not match:
        return None
    library_id = int(match.group(1))
    lines = block.splitlines()
    source_name = lines[0].strip() if lines else ""
    cas_match = _CAS_RN_RE.search("\n".join(lines[:4]))
    source_cas = _cas_from_digits(cas_match.group(1)) if cas_match else None

    Chem, rdBase, Descriptors, rdMolDescriptors, MurckoScaffold = _rdkit()
    with rdBase.BlockLogs():
        mol = Chem.MolFromMolBlock(
            _molblock(block),
            sanitize=True,
            removeHs=False,
            strictParsing=False,
        )
    if mol is None:
        return NistStructureEntry(
            library_id=library_id,
            source_name=source_name,
            source_cas=source_cas,
            canonical_smiles=None,
            inchikey=None,
            molecular_formula=None,
            exact_mass=None,
            murcko_scaffold_smiles=None,
            scaffold_id=None,
            parse_error="rdkit_parse_failed",
        )

    canonical_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    inchikey = Chem.MolToInchiKey(mol)
    formula = rdMolDescriptors.CalcMolFormula(mol)
    exact_mass = float(Descriptors.ExactMolWt(mol))
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    scaffold_smiles = Chem.MolToSmiles(scaffold, canonical=True, isomericSmiles=True)
    return NistStructureEntry(
        library_id=library_id,
        source_name=source_name,
        source_cas=source_cas,
        canonical_smiles=canonical_smiles,
        inchikey=inchikey,
        molecular_formula=formula,
        exact_mass=exact_mass,
        murcko_scaffold_smiles=scaffold_smiles,
        scaffold_id=scaffold_id(scaffold_smiles),
        parse_error=None,
    )


def read_swgdrug_nist_structure_archive(
    path: Path,
) -> tuple[list[NistStructureEntry], str, str]:
    """Read library-linked structures from a SWGDRUG NIST-format ZIP."""
    if not zipfile.is_zipfile(path):
        raise ValueError(f"Not a ZIP archive: {path}")
    with zipfile.ZipFile(path) as archive:
        members = [
            item
            for item in archive.infolist()
            if not item.is_dir() and Path(item.filename).name.upper() == "USRSTRUC.DB"
        ]
        if len(members) != 1:
            raise ValueError(
                "Expected exactly one USRSTRUC.DB member in the NIST archive; "
                f"found {len(members)}"
            )
        member = members[0]
        payload = archive.read(member)

    member_sha256 = hashlib.sha256(payload).hexdigest()
    blocks = re.split(br"\$\$\$\$\r?\n", payload)
    entries: list[NistStructureEntry] = []
    for raw_block in blocks:
        if not raw_block.strip():
            continue
        block = raw_block.decode("latin-1", errors="replace").strip("\x00\r\n ")
        entry = _parse_structure_block(block)
        if entry is not None:
            entries.append(entry)
    if not entries:
        raise ValueError("No library-linked NIST structure records found")
    return entries, member.filename, member_sha256


def _record_library_id(record: RawSpectrum) -> int | None:
    value = record.provenance.source_record_id or ""
    match = re.search(r"(\d+)\s*$", value)
    return int(match.group(1)) if match else None


def _record_names(record: RawSpectrum) -> set[str]:
    values = {
        record.compound_name_raw,
        record.metadata.get("CAS_NAME"),
        record.metadata.get("COMPOUND_LIST_NAME"),
    }
    return {
        normalize_compound_name(value).casefold()
        for value in values
        if value and normalize_compound_name(value)
    }


def _entry_name(entry: NistStructureEntry) -> str:
    return normalize_compound_name(entry.source_name).casefold()


def _record_formula(record: RawSpectrum) -> str | None:
    return record.metadata.get("MOLFORM") or record.metadata.get("FORMULA")


def _record_cas(record: RawSpectrum) -> str | None:
    return normalize_cas(record.metadata.get("CAS_REGISTRY_NO") or record.metadata.get("CAS"))


def _compatible(record: RawSpectrum, entry: NistStructureEntry) -> str | None:
    return structure_formula_relation(_record_formula(record), entry.molecular_formula)


def _candidate_priority(method: str) -> int:
    return {"library_id": 0, "name_formula": 1, "cas_formula": 2}.get(method, 99)


def align_swgdrug_nist_structures(
    records: list[RawSpectrum],
    entries: list[NistStructureEntry],
    *,
    archive_member: str,
    member_sha256: str,
    exact_mass_tolerance: float = 1e-4,
) -> tuple[list[RawSpectrum], NistStructureAlignment, list[dict[str, Any]]]:
    """Attach only uniquely supported NIST structures to SWGDRUG records.

    Mapping priority is current library ID, then exact source name + compatible
    formula, then valid CAS + compatible formula. Multiple distinct structures
    targeting one spectrum are quarantined rather than guessed.
    """
    by_id: dict[int, int] = {}
    by_name: dict[str, list[int]] = defaultdict(list)
    by_cas: dict[str, list[int]] = defaultdict(list)
    for index, record in enumerate(records):
        library_id = _record_library_id(record)
        if library_id is not None:
            by_id[library_id] = index
        for name in _record_names(record):
            by_name[name].append(index)
        cas = _record_cas(record)
        if cas:
            by_cas[cas].append(index)

    proposals: dict[int, list[tuple[NistStructureEntry, str, str]]] = defaultdict(list)
    review: list[dict[str, Any]] = []
    parse_failures = 0

    for entry in entries:
        if entry.parse_error:
            parse_failures += 1
            source_index = by_id.get(entry.library_id)
            review.append(
                {
                    "source_record_id": (
                        records[source_index].provenance.source_record_id or ""
                        if source_index is not None
                        else ""
                    ),
                    "nist_library_id": entry.library_id,
                    "reason": entry.parse_error,
                    "match_method": "library_id" if source_index is not None else "",
                    "source_formula": (
                        _record_formula(records[source_index]) or ""
                        if source_index is not None
                        else ""
                    ),
                    "structure_formula": "",
                }
            )
            continue

        direct_index = by_id.get(entry.library_id)
        if direct_index is not None:
            relation = _compatible(records[direct_index], entry)
            if relation:
                proposals[direct_index].append((entry, "library_id", relation))
                continue

        name_candidates = {
            index
            for index in by_name.get(_entry_name(entry), [])
            if _compatible(records[index], entry)
        }
        if len(name_candidates) == 1:
            index = next(iter(name_candidates))
            relation = _compatible(records[index], entry)
            assert relation is not None
            proposals[index].append((entry, "name_formula", relation))
            continue

        if entry.source_cas:
            cas_candidates = {
                index
                for index in by_cas.get(entry.source_cas, [])
                if _compatible(records[index], entry)
            }
            if len(cas_candidates) == 1:
                index = next(iter(cas_candidates))
                relation = _compatible(records[index], entry)
                assert relation is not None
                proposals[index].append((entry, "cas_formula", relation))
                continue

        if direct_index is not None:
            review.append(
                {
                    "source_record_id": records[direct_index].provenance.source_record_id or "",
                    "nist_library_id": entry.library_id,
                    "reason": "library_id_formula_mismatch_or_ambiguous_fallback",
                    "match_method": "library_id",
                    "source_formula": _record_formula(records[direct_index]) or "",
                    "structure_formula": entry.molecular_formula or "",
                }
            )

    aligned = list(records)
    method_counts: Counter[str] = Counter()
    relation_counts: Counter[str] = Counter()
    matched = 0
    conflicts = 0
    duplicate_collapses = 0
    exact_mass_checked = 0
    exact_mass_warnings = 0

    for index, candidates in proposals.items():
        smiles = {entry.canonical_smiles for entry, _, _ in candidates}
        if len(smiles) != 1:
            conflicts += 1
            for entry, method, relation in candidates:
                review.append(
                    {
                        "source_record_id": records[index].provenance.source_record_id or "",
                        "nist_library_id": entry.library_id,
                        "reason": "multiple_distinct_structure_candidates",
                        "match_method": method,
                        "source_formula": _record_formula(records[index]) or "",
                        "structure_formula": entry.molecular_formula or "",
                        "formula_relation": relation,
                    }
                )
            continue
        if len(candidates) > 1:
            duplicate_collapses += len(candidates) - 1
        entry, method, relation = min(
            candidates,
            key=lambda item: (_candidate_priority(item[1]), item[0].library_id),
        )
        metadata = dict(records[index].metadata)
        metadata.update(
            {
                "NIST_STRUCTURE_CANONICAL_SMILES": entry.canonical_smiles or "",
                "NIST_STRUCTURE_INCHIKEY": entry.inchikey or "",
                "NIST_STRUCTURE_MOLFORM": entry.molecular_formula or "",
                "NIST_STRUCTURE_EXACT_MASS": (
                    f"{entry.exact_mass:.12g}" if entry.exact_mass is not None else ""
                ),
                "NIST_STRUCTURE_SCAFFOLD_SMILES": entry.murcko_scaffold_smiles or "",
                "NIST_STRUCTURE_SCAFFOLD_ID": entry.scaffold_id or "",
                "NIST_STRUCTURE_LIBRARY_ID": str(entry.library_id),
                "NIST_STRUCTURE_MATCH_METHOD": method,
                "NIST_STRUCTURE_FORMULA_RELATION": relation,
                "NIST_STRUCTURE_SOURCE_MEMBER": archive_member,
                "NIST_STRUCTURE_SOURCE_MEMBER_SHA256": member_sha256,
            }
        )
        source_exact = metadata.get("EXACT_MASS_LIST_MONOISOTOPIC_MASS")
        if relation == "exact_formula" and source_exact and entry.exact_mass is not None:
            try:
                delta = abs(float(source_exact) - entry.exact_mass)
            except ValueError:
                pass
            else:
                exact_mass_checked += 1
                metadata["NIST_STRUCTURE_EXACT_MASS_DELTA"] = f"{delta:.12g}"
                if delta > exact_mass_tolerance:
                    exact_mass_warnings += 1
                    metadata["NIST_STRUCTURE_EXACT_MASS_WARNING"] = (
                        f"delta_gt_{exact_mass_tolerance:g}"
                    )
                    review.append(
                        {
                            "source_record_id": (
                                records[index].provenance.source_record_id or ""
                            ),
                            "nist_library_id": entry.library_id,
                            "reason": "exact_mass_delta_warning",
                            "match_method": method,
                            "source_formula": _record_formula(records[index]) or "",
                            "structure_formula": entry.molecular_formula or "",
                            "formula_relation": relation,
                            "source_exact_mass": source_exact,
                            "structure_exact_mass": f"{entry.exact_mass:.12g}",
                            "exact_mass_delta": f"{delta:.12g}",
                        }
                    )
        aligned[index] = replace(records[index], metadata=metadata)
        matched += 1
        method_counts[method] += 1
        relation_counts[relation] += 1

    alignment = NistStructureAlignment(
        linked_structure_blocks=len(entries),
        parsed_structure_blocks=len(entries) - parse_failures,
        parse_failures=parse_failures,
        matched_spectra=matched,
        unmatched_spectra=len(records) - matched,
        conflicting_spectra=conflicts,
        duplicate_candidate_collapses=duplicate_collapses,
        match_method_counts=dict(sorted(method_counts.items())),
        formula_relation_counts=dict(sorted(relation_counts.items())),
        exact_mass_checked_records=exact_mass_checked,
        exact_mass_warning_records=exact_mass_warnings,
    )
    return aligned, alignment, review
