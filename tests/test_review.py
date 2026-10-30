import csv
from pathlib import Path

from opioidbench.chemistry.review import prioritize_taxonomy_review


def test_similarity_review_prioritizes_without_assigning_labels(tmp_path: Path) -> None:
    compounds = tmp_path / "compounds.csv"
    compounds.write_text(
        "compound_id,preferred_name,cas,canonical_smiles,scaffold_id\n"
        "id:opioid,Reference opioid,,c1ccccc1,murcko:benzene\n"
        "id:close,Close analogue,,Cc1ccccc1,murcko:benzene\n"
        "id:far,Dissimilar,,CCCCCC,murcko:acyclic\n",
        encoding="utf-8",
    )
    taxonomy = tmp_path / "taxonomy.csv"
    taxonomy.write_text(
        "compound_id,preferred_name,opioid_status,opioid_family,subfamily,"
        "evidence_source,curation_status,notes\n"
        "id:opioid,Reference opioid,opioid,test,test,test,curated,test\n",
        encoding="utf-8",
    )
    output = tmp_path / "priority.csv"
    counts = prioritize_taxonomy_review(
        compounds,
        taxonomy,
        output_path=output,
        threshold=0.2,
    )
    rows = list(csv.DictReader(output.open()))
    assert counts["prioritized_for_review"] == 1
    assert [row["compound_id"] for row in rows] == ["id:close"]
    assert rows[0]["review_reason"] == "structure_similarity_to_curated_opioid"
    assert "opioid_status" not in rows[0]


def test_similarity_review_excludes_curated_nonopioid(tmp_path: Path) -> None:
    compounds = tmp_path / "compounds.csv"
    compounds.write_text(
        "compound_id,preferred_name,cas,canonical_smiles,scaffold_id\n"
        "id:opioid,Reference opioid,,c1ccccc1,murcko:benzene\n"
        "id:negative,Curated negative,,Cc1ccccc1,murcko:benzene\n",
        encoding="utf-8",
    )
    taxonomy = tmp_path / "taxonomy.csv"
    taxonomy.write_text(
        "compound_id,preferred_name,opioid_status,opioid_family,subfamily,"
        "evidence_source,curation_status,notes\n"
        "id:opioid,Reference opioid,opioid,test,test,test,curated,test\n"
        "id:negative,Curated negative,non-opioid,non-opioid comparator,test,test,curated,test\n",
        encoding="utf-8",
    )
    output = tmp_path / "priority.csv"
    counts = prioritize_taxonomy_review(
        compounds,
        taxonomy,
        output_path=output,
        threshold=0.0,
    )
    assert counts["unresolved_compounds"] == 0
    assert counts["prioritized_for_review"] == 0


def test_similarity_review_excludes_target_excluded_rows(tmp_path: Path) -> None:
    compounds = tmp_path / "compounds.csv"
    compounds.write_text(
        "compound_id,preferred_name,cas,canonical_smiles,scaffold_id\n"
        "id:opioid,Reference opioid,,c1ccccc1,murcko:benzene\n"
        "id:excluded,Target exclusion,,Cc1ccccc1,murcko:benzene\n",
        encoding="utf-8",
    )
    taxonomy = tmp_path / "taxonomy.csv"
    taxonomy.write_text(
        "compound_id,preferred_name,opioid_status,opioid_family,subfamily,"
        "evidence_source,curation_status,notes\n"
        "id:opioid,Reference opioid,opioid,test,test,test,curated,test\n"
        "id:excluded,Target exclusion,,benchmark exclusion,test,test,excluded,test\n",
        encoding="utf-8",
    )
    output = tmp_path / "priority.csv"
    counts = prioritize_taxonomy_review(
        compounds,
        taxonomy,
        output_path=output,
        threshold=0.0,
    )
    assert counts["unresolved_compounds"] == 0
    assert counts["prioritized_for_review"] == 0
