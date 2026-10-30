from pathlib import Path

import pytest

from opioidbench.chemistry.taxonomy import load_taxonomy


def test_empty_taxonomy_header_is_valid(tmp_path: Path) -> None:
    path = tmp_path / "taxonomy.csv"
    path.write_text(
        "compound_id,preferred_name,opioid_status,opioid_family,subfamily,"
        "evidence_source,curation_status,notes\n",
        encoding="utf-8",
    )
    assert load_taxonomy(path) == {}


def test_taxonomy_rejects_uncontrolled_status(tmp_path: Path) -> None:
    path = tmp_path / "taxonomy.csv"
    path.write_text(
        "compound_id,preferred_name,opioid_status,opioid_family,subfamily,"
        "evidence_source,curation_status,notes\n"
        "cas:1-23-4,X,uncertain,,,,guessed,\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="curation_status"):
        load_taxonomy(path)
