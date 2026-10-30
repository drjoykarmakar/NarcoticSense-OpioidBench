import csv
from pathlib import Path

from opioidbench.chemistry.structures import calculate_structure_features
from opioidbench.data.external import (
    _load_supplemental_structures,
    _supplemental_structure_for_row,
)


def test_supplemental_structure_requires_exact_name_and_formula(tmp_path: Path) -> None:
    features = calculate_structure_features("CCO")
    path = tmp_path / "structures.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "preferred_name",
                "aliases",
                "molecular_formula",
                "canonical_smiles",
                "inchikey",
                "pubchem_cid",
                "source",
                "source_url",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "preferred_name": "Ethanol",
                "aliases": "ethyl alcohol",
                "molecular_formula": "C2H6O",
                "canonical_smiles": "CCO",
                "inchikey": features.inchikey,
                "pubchem_cid": "702",
                "source": "PubChem",
                "source_url": "https://example.invalid/702",
            }
        )

    index, count = _load_supplemental_structures(path)
    assert count == 1
    match = _supplemental_structure_for_row(
        {
            "preferred_name": "ethyl alcohol",
            "aliases": "",
            "molecular_formula": "C2H6O",
        },
        index,
    )
    assert match is not None
    assert match.inchikey == features.inchikey

    assert (
        _supplemental_structure_for_row(
            {
                "preferred_name": "ethyl alcohol",
                "aliases": "",
                "molecular_formula": "C2H4O2",
            },
            index,
        )
        is None
    )
    assert (
        _supplemental_structure_for_row(
            {
                "preferred_name": "methanol",
                "aliases": "",
                "molecular_formula": "C2H6O",
            },
            index,
        )
        is None
    )


def test_repository_and_packaged_external_structure_seed_are_identical() -> None:
    root = Path(__file__).parents[1]
    assert (root / "metadata" / "external_structure_evidence.csv").read_bytes() == (
        root / "src" / "opioidbench" / "resources" / "external_structure_evidence.csv"
    ).read_bytes()
