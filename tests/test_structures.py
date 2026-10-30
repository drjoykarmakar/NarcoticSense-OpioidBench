import pytest

pytest.importorskip("rdkit")

from opioidbench.chemistry.structures import calculate_structure_features, morgan_fingerprint


def test_rdkit_structure_features_and_fingerprint() -> None:
    features = calculate_structure_features("CCO")
    assert features.canonical_smiles == "CCO"
    assert len(features.inchikey) == 27
    assert features.molecular_weight > 40
    fp = morgan_fingerprint("CCO", radius=2, n_bits=512, use_chirality=True)
    assert fp.GetNumBits() == 512


def test_invalid_smiles_is_rejected() -> None:
    with pytest.raises(ValueError, match="Invalid SMILES"):
        calculate_structure_features("not a smiles")


def test_scaffold_id_is_stable_and_marks_acyclic_molecules() -> None:
    from opioidbench.chemistry.structures import scaffold_id

    assert scaffold_id("") == "murcko:acyclic"
    assert scaffold_id("c1ccccc1") == scaffold_id("c1ccccc1")
    assert scaffold_id("c1ccccc1").startswith("murcko:")
