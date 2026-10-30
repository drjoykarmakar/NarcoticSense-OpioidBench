import gzip
from pathlib import Path

import pytest

pytest.importorskip("rdkit")

from opioidbench.chemistry.chebi import iter_chebi_sdf, resolve_chebi_compound_table

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

>  <SYNONYMS>
ethyl alcohol

$$$$
"""


def test_chebi_sdf_streams_gzip_and_resolves_exact_cas(tmp_path: Path) -> None:
    sdf = tmp_path / "chebi.sdf.gz"
    with gzip.open(sdf, "wt", encoding="utf-8") as handle:
        handle.write(SDF)
    entries = list(iter_chebi_sdf(sdf))
    assert len(entries) == 1
    assert entries[0].chebi_id == "CHEBI:16236"
    assert entries[0].cas_numbers == ("64-17-5",)

    compounds = tmp_path / "compounds.csv"
    compounds.write_text(
        "compound_id,preferred_name,aliases,cas,molecular_formula,source_monoisotopic_mass,"
        "inchikey,canonical_smiles\n"
        "cas:64-17-5,Ethanol,Ethanol,64-17-5,C2H6O,46.041864812,,\n",
        encoding="utf-8",
    )
    matches, counts = resolve_chebi_compound_table(compounds, chebi_sdf=sdf)
    assert counts["resolved"] == 1
    assert matches[0].chebi_id == "CHEBI:16236"
    assert matches[0].match_method == "cas+name"
    assert matches[0].canonical_smiles == "CCO"


def test_chebi_name_only_match_requires_formula_evidence(tmp_path: Path) -> None:
    sdf = tmp_path / "chebi.sdf"
    sdf.write_text(SDF, encoding="utf-8")
    compounds = tmp_path / "compounds.csv"
    compounds.write_text(
        "compound_id,preferred_name,aliases,cas,molecular_formula,source_monoisotopic_mass,"
        "inchikey,canonical_smiles\n"
        "name:x,Ethanol,Ethanol,,,46.041864812,,\n",
        encoding="utf-8",
    )
    matches, _ = resolve_chebi_compound_table(compounds, chebi_sdf=sdf)
    assert matches[0].status == "conflict"
    assert "name_match_missing_formula_evidence" in (matches[0].error or "")


def test_chebi_cas_match_does_not_override_source_identity_conflict(tmp_path: Path) -> None:
    sdf = tmp_path / "chebi.sdf"
    sdf.write_text(SDF, encoding="utf-8")
    compounds = tmp_path / "compounds.csv"
    compounds.write_text(
        "compound_id,preferred_name,aliases,cas,molecular_formula,source_monoisotopic_mass,"
        "inchikey,canonical_smiles,name_conflict,formula_conflict,exact_mass_conflict,"
        "structure_conflict\n"
        "cas:64-17-5,,,64-17-5,,,,,True,True,True,True\n",
        encoding="utf-8",
    )
    matches, _ = resolve_chebi_compound_table(compounds, chebi_sdf=sdf)
    assert matches[0].status == "conflict"
    assert "source_identity_conflict" in (matches[0].error or "")
