import csv
import json
from pathlib import Path

import pytest

pytest.importorskip("rdkit")

from opioidbench.chemistry.enrichment import enrich_compound_table


def test_enrichment_table_is_offline_and_noninventive_by_default(tmp_path: Path) -> None:
    source = tmp_path / "compounds.csv"
    with source.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["compound_id", "cas"])
        writer.writeheader()
        writer.writerow({"compound_id": "cas:64-17-5", "cas": "64-17-5"})
        writer.writerow({"compound_id": "name:abc", "cas": ""})
    output = tmp_path / "out"
    counts = enrich_compound_table(
        source,
        output_dir=output,
        cache_dir=tmp_path / "cache",
    )
    assert counts == {
        "input_compounds": 2,
        "resolved": 0,
        "missing_or_invalid_cas": 1,
        "cache_miss_offline": 1,
        "not_found": 0,
        "source_formula_mismatch": 0,
        "source_exact_mass_mismatch": 0,
        "errors": 0,
    }
    report = json.loads((output / "enrichment_report.json").read_text(encoding="utf-8"))
    assert report == counts
    descriptors = (output / "calculated_descriptors.csv").read_text(encoding="utf-8").splitlines()
    assert len(descriptors) == 1


def test_enrichment_requires_source_formula_and_exact_mass_consistency(tmp_path: Path) -> None:
    from opioidbench.chemistry.enrichment import MetadataCache

    source = tmp_path / "compounds.csv"
    fieldnames = [
        "compound_id",
        "cas",
        "molecular_formula",
        "source_monoisotopic_mass",
    ]
    with source.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(
            {
                "compound_id": "cas:64-17-5",
                "cas": "64-17-5",
                "molecular_formula": "C2H6O",
                "source_monoisotopic_mass": "46.041864812",
            }
        )

    cache = MetadataCache(tmp_path / "cache")
    cache.put(
        "pubchem",
        "cas",
        "64-17-5",
        {
            "status": "resolved",
            "record": {
                "query_type": "cas",
                "query": "64-17-5",
                "source": "PubChem PUG REST",
                "source_url": "https://example.invalid/pubchem",
                "retrieval_date": "2026-09-29",
                "pubchem_cid": "702",
                "smiles": "CCO",
                "canonical_smiles": "CCO",
                "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                "molecular_formula": "C2H6O",
                "molecular_weight": 46.07,
                "exact_mass": 46.041864812,
            },
        },
    )
    output = tmp_path / "out"
    counts = enrich_compound_table(source, output_dir=output, cache_dir=tmp_path / "cache")
    assert counts["resolved"] == 1
    assert counts["source_formula_mismatch"] == 0
    assert counts["source_exact_mass_mismatch"] == 0
    assert len((output / "calculated_descriptors.csv").read_text().splitlines()) == 2


def test_enrichment_quarantines_formula_mismatch(tmp_path: Path) -> None:
    from opioidbench.chemistry.enrichment import MetadataCache

    source = tmp_path / "compounds.csv"
    with source.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["compound_id", "cas", "molecular_formula", "source_monoisotopic_mass"],
        )
        writer.writeheader()
        writer.writerow(
            {
                "compound_id": "cas:64-17-5",
                "cas": "64-17-5",
                "molecular_formula": "C3H8O",
                "source_monoisotopic_mass": "60.057514876",
            }
        )
    cache = MetadataCache(tmp_path / "cache")
    cache.put(
        "pubchem",
        "cas",
        "64-17-5",
        {
            "status": "resolved",
            "record": {
                "query_type": "cas",
                "query": "64-17-5",
                "source": "PubChem PUG REST",
                "source_url": "https://example.invalid/pubchem",
                "retrieval_date": "2026-09-29",
                "pubchem_cid": "702",
                "smiles": "CCO",
                "canonical_smiles": "CCO",
                "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                "molecular_formula": "C2H6O",
                "molecular_weight": 46.07,
                "exact_mass": 46.041864812,
            },
        },
    )
    output = tmp_path / "out"
    counts = enrich_compound_table(source, output_dir=output, cache_dir=tmp_path / "cache")
    assert counts["resolved"] == 0
    assert counts["source_formula_mismatch"] == 1
    assert len((output / "calculated_descriptors.csv").read_text().splitlines()) == 1
