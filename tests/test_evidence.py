import csv
from pathlib import Path

from opioidbench.chemistry.evidence import curate_taxonomy_from_evidence
from opioidbench.chemistry.taxonomy import load_taxonomy


def test_evidence_curation_marks_only_exact_matches_and_keeps_others_pending(
    tmp_path: Path,
) -> None:
    compounds = tmp_path / "compounds.csv"
    compounds.write_text(
        "compound_id,preferred_name,aliases,cas,molecular_formula,inchikey,canonical_smiles\n"
        "cas:57-27-2,Morphine,Morphine,57-27-2,C17H19NO3,,\n"
        "name:x,Fentanylish,Fentanylish,,C1H1,,\n",
        encoding="utf-8",
    )
    evidence = tmp_path / "evidence.csv"
    fields = [
        "evidence_id",
        "preferred_name",
        "aliases",
        "cas",
        "opioid_status",
        "opioid_family",
        "subfamily",
        "evidence_source",
        "evidence_url",
        "evidence_date",
        "identity_source",
        "identity_url",
        "notes",
    ]
    with evidence.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerow(
            {
                "evidence_id": "x1",
                "preferred_name": "morphine",
                "aliases": "",
                "cas": "",
                "opioid_status": "opioid",
                "opioid_family": "natural opium alkaloids",
                "subfamily": "N02AA",
                "evidence_source": "WHO ATC",
                "evidence_url": "https://example.invalid",
                "evidence_date": "2026-01-20",
                "identity_source": "",
                "identity_url": "",
                "notes": "",
            }
        )
    taxonomy = tmp_path / "taxonomy.csv"
    review = tmp_path / "review.csv"
    counts = curate_taxonomy_from_evidence(
        compounds,
        evidence,
        taxonomy_output=taxonomy,
        review_output=review,
    )
    assert counts["curated_opioid"] == 1
    assert counts["pending"] == 1
    loaded = load_taxonomy(taxonomy)
    assert loaded["cas:57-27-2"].opioid_status == "opioid"
    assert "name:x" not in loaded
    assert "Fentanylish" in review.read_text(encoding="utf-8")


def test_repository_and_packaged_evidence_seed_are_identical() -> None:
    root = Path(__file__).parents[1]
    assert (root / "metadata" / "opioid_evidence.csv").read_bytes() == (
        root / "src" / "opioidbench" / "resources" / "opioid_evidence.csv"
    ).read_bytes()


def test_evidence_curation_supports_explicit_nonopioid_comparators(tmp_path: Path) -> None:
    compounds = tmp_path / "compounds.csv"
    compounds.write_text(
        "compound_id,preferred_name,aliases,cas,molecular_formula,inchikey,canonical_smiles\n"
        "cas:439-14-5,Diazepam,Diazepam,439-14-5,C16H13ClN2O,,\n",
        encoding="utf-8",
    )
    evidence = tmp_path / "evidence.csv"
    fields = [
        "evidence_id",
        "preferred_name",
        "aliases",
        "cas",
        "opioid_status",
        "opioid_family",
        "subfamily",
        "evidence_source",
        "evidence_url",
        "evidence_date",
        "identity_source",
        "identity_url",
        "notes",
    ]
    with evidence.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerow(
            {
                "evidence_id": "n1",
                "preferred_name": "diazepam",
                "opioid_status": "non-opioid",
                "opioid_family": "non-opioid comparator",
                "subfamily": "benzodiazepine anxiolytic (N05BA)",
                "evidence_source": "WHO ATC",
                "evidence_url": "https://example.invalid",
                "evidence_date": "2026-01-20",
            }
        )
    taxonomy = tmp_path / "taxonomy.csv"
    review = tmp_path / "review.csv"
    counts = curate_taxonomy_from_evidence(
        compounds,
        evidence,
        taxonomy_output=taxonomy,
        review_output=review,
    )
    assert counts["curated_nonopioid"] == 1
    loaded = load_taxonomy(taxonomy)
    assert loaded["cas:439-14-5"].opioid_status == "non-opioid"


def test_conflicting_exact_evidence_goes_to_uncertain_review(tmp_path: Path) -> None:
    compounds = tmp_path / "compounds.csv"
    compounds.write_text(
        "compound_id,preferred_name,aliases,cas,molecular_formula,inchikey,canonical_smiles\n"
        "cas:1-11-1,Example,Example,,C1H1,,\n",
        encoding="utf-8",
    )
    evidence = tmp_path / "evidence.csv"
    fields = [
        "evidence_id",
        "preferred_name",
        "aliases",
        "cas",
        "opioid_status",
        "opioid_family",
        "subfamily",
        "evidence_source",
        "evidence_url",
        "evidence_date",
        "identity_source",
        "identity_url",
        "notes",
    ]
    with evidence.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerow(
            {
                "evidence_id": "p1",
                "preferred_name": "Example",
                "opioid_status": "opioid",
                "opioid_family": "test",
                "evidence_source": "source A",
            }
        )
        writer.writerow(
            {
                "evidence_id": "n1",
                "preferred_name": "Example",
                "opioid_status": "non-opioid",
                "opioid_family": "non-opioid comparator",
                "evidence_source": "source B",
            }
        )
    taxonomy = tmp_path / "taxonomy.csv"
    review = tmp_path / "review.csv"
    counts = curate_taxonomy_from_evidence(
        compounds,
        evidence,
        taxonomy_output=taxonomy,
        review_output=review,
    )
    assert counts["uncertain_conflict"] == 1
    assert load_taxonomy(taxonomy) == {}
    assert "conflicting_authoritative_evidence" in review.read_text(encoding="utf-8")


def test_evidence_curation_excludes_analytical_variants_and_target_exclusions(
    tmp_path: Path,
) -> None:
    compounds = tmp_path / "compounds.csv"
    compounds.write_text(
        "compound_id,preferred_name,aliases,cas,molecular_formula,inchikey,canonical_smiles\n"
        "id:d3,Oxycodone-d3,,145225-02-1,C18H18D3NO4,,[2H]C\n"
        "id:tms,Codeine TMS Derivative,,74367-14-9,C21H29NO3Si,,CC\n"
        "id:ant,Dextrorphan,,125-73-5,C17H23NO,,CC\n",
        encoding="utf-8",
    )
    evidence = tmp_path / "evidence.csv"
    evidence.write_text(
        ",".join(
            [
                "evidence_id",
                "preferred_name",
                "aliases",
                "cas",
                "opioid_status",
                "opioid_family",
                "subfamily",
                "evidence_source",
                "evidence_url",
                "evidence_date",
                "identity_source",
                "identity_url",
                "notes",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    exclusions = tmp_path / "exclusions.csv"
    exclusions.write_text(
        "exclusion_id,preferred_name,aliases,cas,reason,evidence_source,evidence_url,"
        "evidence_date,notes\n"
        "x1,Dextrorphan,,125-73-5,opioid_related_nonagonist_outside_target,test,"
        "https://example.invalid,2026,test\n",
        encoding="utf-8",
    )
    taxonomy = tmp_path / "taxonomy.csv"
    review = tmp_path / "review.csv"
    counts = curate_taxonomy_from_evidence(
        compounds,
        evidence,
        taxonomy_output=taxonomy,
        review_output=review,
        exclusions_path=exclusions,
    )
    assert counts["excluded_analytical_variant"] == 2
    assert counts["excluded_target"] == 1
    assert counts["pending"] == 0
    loaded = load_taxonomy(taxonomy)
    assert set(loaded) == {"id:d3", "id:tms", "id:ant"}
    assert all(entry.curation_status == "excluded" for entry in loaded.values())
    assert loaded["id:ant"].opioid_status == ""
    assert review.read_text(encoding="utf-8").count("\n") == 1


def test_repository_and_packaged_exclusion_seed_are_identical() -> None:
    root = Path(__file__).parents[1]
    assert (root / "metadata" / "taxonomy_exclusions.csv").read_bytes() == (
        root / "src" / "opioidbench" / "resources" / "taxonomy_exclusions.csv"
    ).read_bytes()


def test_repository_nbenzyl_fentanyl_edge_cases_are_target_excluded() -> None:
    from opioidbench.chemistry.exclusions import (
        load_taxonomy_exclusions,
        match_taxonomy_exclusion,
    )

    root = Path(__file__).parents[1]
    exclusions = load_taxonomy_exclusions(root / "metadata" / "taxonomy_exclusions.csv")
    cases = [
        ("Benzyl Carfentanil", "61085-72-1"),
        ("Benzyl Acrylfentanyl", "1956368-29-8"),
        ("N-benzyl acetyl fentanyl", "1237-52-1"),
    ]
    for name, cas in cases:
        match = match_taxonomy_exclusion(
            {"preferred_name": name, "aliases": "", "cas": cas}, exclusions
        )
        assert match is not None
        assert match.reason == "fentanyl_related_nbenzyl_pharmacology_unresolved"


def test_repository_review_frontier_evidence_has_exact_aliases() -> None:
    from opioidbench.chemistry.evidence import load_opioid_evidence

    root = Path(__file__).parents[1]
    evidence = load_opioid_evidence(root / "metadata" / "opioid_evidence.csv")
    by_id = {row.evidence_id: row for row in evidence}
    assert "6-Monoacetylmorphine" in by_id["review51-frontier-6mam"].aliases
    assert "Fentanyl ortho tolyl acetyl analog" in (
        by_id["review51-frontier-ortho-methyl-acetyl-fentanyl"].aliases
    )
    assert by_id["review51-frontier-ortho-fluoro-acrylfentanyl"].opioid_status == "opioid"
