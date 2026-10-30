import pytest

from opioidbench.chemistry.identifiers import (
    cas_checksum_valid,
    conservative_compound_id,
    inchikey_valid,
    molecular_formula_plausible,
    normalize_cas,
    normalize_compound_name,
)


def test_cas_normalization_and_placeholders() -> None:
    assert cas_checksum_valid("000481-29-8")
    assert normalize_cas("000481-29-8") == "481-29-8"
    assert normalize_cas("000000-00-0") is None
    assert not cas_checksum_valid("185888-30-5")


def test_name_normalization_is_conservative() -> None:
    assert normalize_compound_name("  Alpha   Beta  ") == "Alpha Beta"
    assert conservative_compound_id(cas=None, name="Alpha Beta").startswith("name:")
    assert conservative_compound_id(cas="000481-29-8", name="Anything") == "cas:481-29-8"
    with pytest.raises(ValueError):
        conservative_compound_id(cas=None, name=None)


def test_identifier_and_formula_syntax_checks() -> None:
    assert inchikey_valid("LFQSCWFLJHTTHZ-UHFFFAOYSA-N")
    assert not inchikey_valid("not-an-inchikey")
    assert molecular_formula_plausible("C19H30O2")
    assert not molecular_formula_plausible("C19 H30 O2")


def test_enfsi_zero_cas_is_treated_as_missing_placeholder() -> None:
    from opioidbench.chemistry.identifiers import CAS_PLACEHOLDERS, normalize_cas

    assert "0" in CAS_PLACEHOLDERS
    assert normalize_cas("0") is None


def test_structure_formula_relation_preserves_deuterated_parent_distinction() -> None:
    from opioidbench.chemistry.identifiers import (
        formula_atom_counts,
        structure_formula_relation,
    )

    assert formula_atom_counts("C22H26NO3+") == {"C": 22, "H": 26, "N": 1, "O": 3}
    assert structure_formula_relation("C2H6O", "C2H6O") == "exact_formula"
    assert structure_formula_relation("C2HD5O", "C2H6O") == "isotopologue_parent"
    assert structure_formula_relation("C2H6O", "C3H8O") is None
