import csv
from pathlib import Path

import pytest

pytest.importorskip("rdkit")

from opioidbench.chemistry.resolution import resolve_compound_sources


SDF = """ethanol
  ChEBI

>  <ChEBI ID>
CHEBI:16236

>  <ChEBI NAME>
ethanol

>  <FORMULA>
C2H6O

>  <MONOISOTOPIC MASS>
46.041864812

>  <SMILES>
CCO

>  <INCHIKEY>
LFQSCWFLJHTTHZ-UHFFFAOYSA-N

>  <CAS Registry Numbers>
64-17-5

$$$$
"""


def test_resolution_is_chebi_first_and_fills_missing_structure(tmp_path: Path) -> None:
    source = tmp_path / "compounds.csv"
    source.write_text(
        "compound_id,preferred_name,aliases,cas,molecular_formula,source_monoisotopic_mass,"
        "inchikey,canonical_smiles,scaffold_id,murcko_scaffold_smiles,structure_source,"
        "structure_formula_relation,pubchem_cid\n"
        "cas:64-17-5,Ethanol,Ethanol,64-17-5,C2H6O,46.041864812,,,,,,,\n",
        encoding="utf-8",
    )
    chebi = tmp_path / "chebi.sdf"
    chebi.write_text(SDF, encoding="utf-8")
    output = tmp_path / "out"
    counts = resolve_compound_sources(
        source,
        output_dir=output,
        chebi_sdf=chebi,
        pubchem_cache_dir=tmp_path / "cache",
    )
    assert counts["chebi_resolved"] == 1
    assert counts["chebi_structure_filled"] == 1
    assert counts["final_structures"] == 1
    assert counts["final_structureless"] == 0
    assert counts["final_unique_scaffold_groups"] == 1
    assert counts["new_scaffold_groups"] == 1
    with (output / "resolved_compounds.csv").open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle))
    assert row["chebi_id"] == "CHEBI:16236"
    assert row["canonical_smiles"] == "CCO"
    assert row["external_resolution_status"] == "chebi_structure_filled"
    assert row["pubchem_cid"] == ""


def test_pubchem_validation_rejects_conflicted_cas_without_existing_inchikey() -> None:
    from opioidbench.chemistry.enrichment import EnrichmentRecord
    from opioidbench.chemistry.resolution import _validate_pubchem_record

    row = {
        "cas": "64-17-5",
        "name_conflict": "True",
        "formula_conflict": "True",
        "structure_conflict": "True",
    }
    record = EnrichmentRecord(
        query_type="cas",
        query="64-17-5",
        source="PubChem PUG REST",
        source_url="https://example.invalid",
        retrieval_date="2026-09-29",
        pubchem_cid="702",
        smiles="CCO",
        canonical_smiles="CCO",
        inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        molecular_formula="C2H6O",
        molecular_weight=46.07,
        exact_mass=46.041864812,
    )
    ok, _, error = _validate_pubchem_record(row, record)
    assert ok is False
    assert error == "source_identity_conflict"


def test_resolution_preserves_existing_source_scaffold(tmp_path: Path) -> None:
    source = tmp_path / "compounds.csv"
    source.write_text(
        "compound_id,preferred_name,aliases,cas,molecular_formula,source_monoisotopic_mass,"
        "inchikey,canonical_smiles,scaffold_id,murcko_scaffold_smiles,structure_source,"
        "structure_formula_relation,pubchem_cid\n"
        "cas:64-17-5,Ethanol,Ethanol,64-17-5,C2H6O,46.041864812,"
        "LFQSCWFLJHTTHZ-UHFFFAOYSA-N,CCO,murcko:source-stable,,source-linked,"
        "exact_formula,\n",
        encoding="utf-8",
    )
    chebi = tmp_path / "chebi.sdf"
    chebi.write_text(SDF, encoding="utf-8")
    output = tmp_path / "out"
    counts = resolve_compound_sources(
        source,
        output_dir=output,
        chebi_sdf=chebi,
        pubchem_cache_dir=tmp_path / "cache",
    )
    assert counts["chebi_resolved"] == 1
    with (output / "resolved_compounds.csv").open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle))
    assert row["external_resolution_status"] == "chebi_identity_resolved_existing_structure"
    assert row["scaffold_id"] == "murcko:source-stable"
    assert row["structure_source"] == "source-linked"
