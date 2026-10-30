import csv
import json
from pathlib import Path

from opioidbench.experiments.splits import generate_cv_splits


def _write_fixture(tmp_path: Path):
    compounds = tmp_path / "compounds.csv"
    with compounds.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["compound_id", "scaffold_id", "structure_conflict"])
        w.writeheader()
        for i in range(12):
            w.writerow(
                {
                    "compound_id": f"c{i}",
                    "scaffold_id": f"s{i // 2}",
                    "structure_conflict": "False",
                }
            )
    taxonomy = tmp_path / "taxonomy.csv"
    fields = [
        "compound_id",
        "preferred_name",
        "opioid_status",
        "opioid_family",
        "subfamily",
        "evidence_source",
        "curation_status",
        "notes",
    ]
    with taxonomy.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for i in range(12):
            w.writerow(
                {
                    "compound_id": f"c{i}",
                    "preferred_name": f"n{i}",
                    "opioid_status": "opioid" if i % 2 else "non-opioid",
                    "opioid_family": "",
                    "subfamily": "",
                    "evidence_source": "x",
                    "curation_status": "curated",
                    "notes": "",
                }
            )
    dataset = tmp_path / "records.jsonl"
    with dataset.open("w", encoding="utf-8") as f:
        for i in range(12):
            f.write(
                json.dumps(
                    {
                        "spectrum_id": f"sp{i}",
                        "compound_id": f"c{i}",
                        "peaks": [{"mz": 50 + i, "intensity": 100.0}],
                    }
                )
                + "\n"
            )
    return dataset, compounds, taxonomy


def test_scaffold_cv_has_zero_scaffold_and_compound_leakage(tmp_path):
    dataset, compounds, taxonomy = _write_fixture(tmp_path)
    report = generate_cv_splits(
        dataset,
        compounds,
        taxonomy,
        output_dir=tmp_path / "out",
        protocol="scaffold",
        n_splits=3,
        seed=7,
    )
    assert report["input"]["included_spectra"] == 12
    assert all(fold["compound_overlap"] == 0 for fold in report["folds"])
    assert all(fold["scaffold_overlap"] == 0 for fold in report["folds"])


def test_analogue_components_union_scaffolds_and_positive_lineages() -> None:
    from opioidbench.experiments.splits import BenchmarkSample, _analogue_component_groups

    samples = [
        BenchmarkSample(
            spectrum_id="sp0",
            compound_id="c0",
            label=1,
            label_name="opioid",
            scaffold_id="s0",
            structure_conflict=False,
            peaks=((50.0, 100.0),),
            analogue_group="fentanyl_anilidopiperidine",
        ),
        BenchmarkSample(
            spectrum_id="sp1",
            compound_id="c1",
            label=1,
            label_name="opioid",
            scaffold_id="s1",
            structure_conflict=False,
            peaks=((51.0, 100.0),),
            analogue_group="fentanyl_anilidopiperidine",
        ),
        BenchmarkSample(
            spectrum_id="sp2",
            compound_id="c2",
            label=0,
            label_name="non-opioid",
            scaffold_id="s1",
            structure_conflict=False,
            peaks=((52.0, 100.0),),
        ),
        BenchmarkSample(
            spectrum_id="sp3",
            compound_id="c3",
            label=0,
            label_name="non-opioid",
            scaffold_id="s2",
            structure_conflict=False,
            peaks=((53.0, 100.0),),
        ),
    ]
    groups = _analogue_component_groups(samples)
    assert groups[0] == groups[1] == groups[2]
    assert groups[3] != groups[0]


def test_split_loader_ignores_excluded_taxonomy_rows(tmp_path: Path) -> None:
    dataset, compounds, taxonomy = _write_fixture(tmp_path)
    rows = list(csv.DictReader(taxonomy.open(newline="", encoding="utf-8")))
    rows[0]["opioid_status"] = ""
    rows[0]["curation_status"] = "excluded"
    with taxonomy.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    report = generate_cv_splits(
        dataset,
        compounds,
        taxonomy,
        output_dir=tmp_path / "out-excluded",
        protocol="compound",
        n_splits=3,
        seed=7,
    )
    assert report["input"]["included_spectra"] == 11
    assignments = list(
        csv.DictReader((tmp_path / "out-excluded" / "compound_cv3_assignments.csv").open())
    )
    assert all(row["compound_id"] != "c0" for row in assignments)
