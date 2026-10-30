from __future__ import annotations

import zipfile
from pathlib import Path

import pytest

pytest.importorskip("rdkit")

from rdkit import Chem

from opioidbench.data.nist import (
    align_swgdrug_nist_structures,
    read_swgdrug_nist_structure_archive,
)
from opioidbench.data.parsers import parse_jcamp_text
from opioidbench.data.validation import normalize_spectrum, validate_dataset


def _molblock(name: str, library_id: int, *, old_variant: bool = False) -> str:
    mol = Chem.MolFromSmiles("CCO")
    assert mol is not None
    lines = Chem.MolToMolBlock(mol).splitlines()
    lines[0] = name
    lines[2] = f"Library ID = {library_id}"
    if old_variant:
        lines[3] = lines[3].replace("V2000", "").rstrip()
        lines = [line for line in lines if line.strip() != "M  END"]
    return "\r\n".join(lines) + "\r\n$$$$\r\n"


def test_nist_archive_parses_old_and_v2000_records(tmp_path: Path) -> None:
    archive = tmp_path / "nist.zip"
    payload = (_molblock("ETHANOL", 10) + _molblock("ETHANOL-D5", 11, old_variant=True)).encode(
        "latin-1"
    )
    with zipfile.ZipFile(archive, "w") as handle:
        handle.writestr("SWGDRUG 3.14/USRSTRUC.DB", payload)

    entries, member, digest = read_swgdrug_nist_structure_archive(archive)
    assert member.endswith("USRSTRUC.DB")
    assert len(digest) == 64
    assert len(entries) == 2
    assert all(entry.parse_error is None for entry in entries)
    assert {entry.molecular_formula for entry in entries} == {"C2H6O"}


def test_nist_alignment_accepts_exact_and_unlabeled_isotopologue_parent(
    tmp_path: Path,
) -> None:
    archive = tmp_path / "nist.zip"
    payload = (_molblock("ETHANOL", 10) + _molblock("ETHANOL-D5", 11)).encode("latin-1")
    with zipfile.ZipFile(archive, "w") as handle:
        handle.writestr("SWGDRUG 3.14/USRSTRUC.DB", payload)
    entries, member, digest = read_swgdrug_nist_structure_archive(archive)

    raw = parse_jcamp_text(
        """
##TITLE=Library Entry 10
##CAS NAME=Ethanol
##MOLFORM=C2H6O
##XYDATA=(XY..XY)
31,100 45,50
##TITLE=Library Entry 11
##CAS NAME=Ethanol-d5
##MOLFORM=C2HD5O
##XYDATA=(XY..XY)
31,100 45,50
""",
        source="swgdrug_eims",
        source_version="3.14",
    )
    aligned, report, review = align_swgdrug_nist_structures(
        raw,
        entries,
        archive_member=member,
        member_sha256=digest,
    )
    assert report.matched_spectra == 2
    assert report.formula_relation_counts == {
        "exact_formula": 1,
        "isotopologue_parent": 1,
    }
    assert not review

    validated, _ = validate_dataset(aligned)
    normalized = [normalize_spectrum(item) for item in validated]
    assert all(item.canonical_smiles == "CCO" for item in normalized)
    assert all(item.inchikey for item in normalized)
    assert normalized[0].structure_formula_relation == "exact_formula"
    assert normalized[1].structure_formula_relation == "isotopologue_parent"
    assert normalized[0].scaffold_id == "murcko:acyclic"
