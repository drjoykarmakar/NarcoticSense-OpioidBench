"""Structure-prioritized review queues for unresolved taxonomy records."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

from opioidbench.chemistry.structures import ChemistryDependencyError, morgan_fingerprint
from opioidbench.chemistry.taxonomy import load_taxonomy

REVIEW_COLUMNS = (
    "compound_id",
    "preferred_name",
    "cas",
    "max_opioid_tanimoto",
    "nearest_opioid_compound_id",
    "nearest_opioid_name",
    "scaffold_id",
    "same_scaffold_as_nearest",
    "review_reason",
)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def prioritize_taxonomy_review(
    compounds_path: Path,
    taxonomy_path: Path,
    *,
    output_path: Path,
    threshold: float = 0.75,
) -> dict[str, int | float | str]:
    """Rank unresolved structured compounds by similarity to curated opioids.

    Similarity is used only to prioritize human review. It never assigns an
    opioid or non-opioid label.
    """
    if not 0.0 <= threshold <= 1.0:
        raise ValueError("threshold must be between 0 and 1")

    try:
        from rdkit import DataStructs
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise ChemistryDependencyError(
            "RDKit is required for taxonomy review prioritization; install opioidbench[chem]."
        ) from exc

    with compounds_path.open(newline="", encoding="utf-8") as handle:
        compounds = list(csv.DictReader(handle))
    taxonomy = load_taxonomy(taxonomy_path)
    compounds_by_id = {
        (row.get("compound_id") or "").strip(): row
        for row in compounds
        if (row.get("compound_id") or "").strip()
    }

    opioid_references: list[tuple[str, str, str, object]] = []
    curated_opioid = 0
    curated_nonopioid = 0
    for compound_id, entry in taxonomy.items():
        if entry.curation_status != "curated":
            continue
        if entry.opioid_status == "opioid":
            curated_opioid += 1
            row = compounds_by_id.get(compound_id)
            smiles = (row or {}).get("canonical_smiles", "") if row else ""
            if smiles:
                opioid_references.append(
                    (
                        compound_id,
                        entry.preferred_name,
                        (row.get("scaffold_id") or "").strip(),
                        morgan_fingerprint(smiles),
                    )
                )
        elif entry.opioid_status == "non-opioid":
            curated_nonopioid += 1

    candidates: list[dict[str, str]] = []
    unresolved_compounds = 0
    unresolved_structured = 0
    for row in compounds:
        compound_id = (row.get("compound_id") or "").strip()
        entry = taxonomy.get(compound_id)
        if entry and entry.curation_status in {"curated", "excluded"}:
            continue
        unresolved_compounds += 1
        smiles = (row.get("canonical_smiles") or "").strip()
        if not smiles:
            continue
        unresolved_structured += 1
        if not opioid_references:
            continue
        fingerprint = morgan_fingerprint(smiles)
        best = max(
            (
                (
                    float(DataStructs.TanimotoSimilarity(fingerprint, ref_fp)),
                    ref_id,
                    ref_name,
                    ref_scaffold,
                )
                for ref_id, ref_name, ref_scaffold, ref_fp in opioid_references
            ),
            key=lambda item: item[0],
        )
        similarity, ref_id, ref_name, ref_scaffold = best
        if similarity < threshold:
            continue
        scaffold = (row.get("scaffold_id") or "").strip()
        candidates.append(
            {
                "compound_id": compound_id,
                "preferred_name": (row.get("preferred_name") or "").strip(),
                "cas": (row.get("cas") or "").strip(),
                "max_opioid_tanimoto": f"{similarity:.6f}",
                "nearest_opioid_compound_id": ref_id,
                "nearest_opioid_name": ref_name,
                "scaffold_id": scaffold,
                "same_scaffold_as_nearest": (
                    "yes" if scaffold and ref_scaffold and scaffold == ref_scaffold else "no"
                ),
                "review_reason": "structure_similarity_to_curated_opioid",
            }
        )

    candidates.sort(
        key=lambda row: (
            -float(row["max_opioid_tanimoto"]),
            row["preferred_name"].casefold(),
            row["compound_id"],
        )
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(REVIEW_COLUMNS), lineterminator="\n")
        writer.writeheader()
        writer.writerows(candidates)

    return {
        "input_compounds": len(compounds),
        "input_compounds_sha256": _sha256_file(compounds_path),
        "taxonomy_sha256": _sha256_file(taxonomy_path),
        "curated_opioid": curated_opioid,
        "curated_nonopioid": curated_nonopioid,
        "structured_opioid_references": len(opioid_references),
        "unresolved_compounds": unresolved_compounds,
        "unresolved_structured": unresolved_structured,
        "similarity_threshold": threshold,
        "prioritized_for_review": len(candidates),
    }
