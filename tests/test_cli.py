import json
import zipfile
from pathlib import Path

from opioidbench.cli import main

ROOT = Path(__file__).parents[1]


def test_demo_writes_qc_normalized_and_manifest(tmp_path: Path) -> None:
    output = tmp_path / "demo"
    code = main(
        [
            "demo",
            "--fixture",
            str(ROOT / "tests" / "fixtures" / "synthetic_eims.jdx"),
            "--output",
            str(output),
            "--registry",
            str(ROOT / "metadata" / "source_registry.yaml"),
        ]
    )
    assert code == 0
    qc = json.loads((output / "data_qc.json").read_text(encoding="utf-8"))
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    rows = (output / "normalized_records.jsonl").read_text(encoding="utf-8").splitlines()
    assert qc["number_spectra"] == 2
    assert qc["number_invalid"] == 0
    assert manifest["scientific_result"] is False
    assert len(rows) == 2


def test_download_command_fails_safely_for_unresolved_swgdrug_policy(
    tmp_path: Path,
) -> None:
    code = main(
        [
            "data",
            "download",
            "--source",
            "swgdrug_eims",
            "--destination",
            str(tmp_path / "swgdrug.dat"),
            "--registry",
            str(ROOT / "metadata" / "source_registry.yaml"),
        ]
    )
    assert code == 2
    assert not (tmp_path / "swgdrug.dat").exists()


def test_import_command_writes_provenance_sidecar(tmp_path: Path) -> None:
    destination = tmp_path / "raw"
    fixture = ROOT / "tests" / "fixtures" / "synthetic_eims.jdx"
    code = main(
        [
            "data",
            "import",
            "--source",
            "synthetic_fixture",
            "--input",
            str(fixture),
            "--destination",
            str(destination),
            "--retrieval-date",
            "2026-09-01",
            "--registry",
            str(ROOT / "metadata" / "source_registry.yaml"),
        ]
    )
    assert code == 0
    imported = destination / fixture.name
    sidecar = imported.with_suffix(imported.suffix + ".provenance.json")
    payload = json.loads(sidecar.read_text(encoding="utf-8"))
    assert payload["source"] == "synthetic_fixture"
    assert payload["retrieval_date"] == "2026-09-01"
    assert len(payload["sha256"]) == 64


def test_import_command_extracts_single_jcamp_zip_with_archive_provenance(
    tmp_path: Path,
) -> None:
    fixture = ROOT / "tests" / "fixtures" / "synthetic_swgdrug_layout.hpj"
    archive_path = tmp_path / "source.zip"
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.write(fixture, arcname="Synthetic Library.HPJ")
    destination = tmp_path / "raw"
    code = main(
        [
            "data",
            "import",
            "--source",
            "synthetic_fixture",
            "--input",
            str(archive_path),
            "--destination",
            str(destination),
            "--registry",
            str(ROOT / "metadata" / "source_registry.yaml"),
        ]
    )
    assert code == 0
    imported = destination / "Synthetic Library.HPJ"
    payload = json.loads(
        imported.with_suffix(imported.suffix + ".provenance.json").read_text(encoding="utf-8")
    )
    assert payload["source_archive"] == "source.zip"
    assert payload["source_archive_member"] == "Synthetic Library.HPJ"
    assert len(payload["source_archive_sha256"]) == 64
    assert len(payload["sha256"]) == 64


def test_enfsi_overlap_cli_writes_aggregate_and_private_outputs(tmp_path: Path) -> None:
    swg = tmp_path / "swg.jdx"
    enf = tmp_path / "enf.jdx"
    swg.write_text(
        "##TITLE=Library Entry 1\n##CAS NAME=A\n##NPOINTS=2\n"
        "##XYDATA=(XY..XY)\n50,1000\n75,500\n##END=\n",
        encoding="utf-8",
    )
    enf.write_text(
        "##TITLE=Library Entry 1\n##CAS NAME=A\n##NPOINTS=2\n"
        "##XYDATA=(XY..XY)\n50,100\n75,50\n##END=\n",
        encoding="utf-8",
    )
    aggregate = tmp_path / "aggregate.json"
    private = tmp_path / "private.jsonl"
    code = main(
        [
            "data",
            "audit-enfsi-overlap",
            "--enfsi",
            str(enf),
            "--swgdrug",
            str(swg),
            "--output",
            str(aggregate),
            "--private-output",
            str(private),
            "--registry",
            str(ROOT / "metadata" / "source_registry.yaml"),
        ]
    )
    assert code == 0
    payload = json.loads(aggregate.read_text(encoding="utf-8"))
    assert payload["near_identical_best_match_records"] == 1
    assert payload["spectrum_disjoint_candidate_records"] == 0
    assert "CAS NAME" not in private.read_text(encoding="utf-8")
