import json
import zipfile
from pathlib import Path

from openpyxl import Workbook

from opioidbench.data.enfsi import audit_enfsi_entry_list, read_enfsi_entry_list
from opioidbench.data.overlap import (
    audit_enfsi_spectral_overlap,
    read_jcamp_source,
    write_private_match_manifest,
)
from opioidbench.data.parsers import parse_jcamp_text


def _write_enfsi(path: Path) -> None:
    workbook = Workbook()
    sheet = workbook.active
    sheet.append(["#", "Name", "Formula", "MW"])
    sheet.append([1, "BZP; 1-benzylpiperazine", "C11H16N2", 176])
    sheet.append([2, "Another Name", "C2H6O", 46])
    sheet.append([3, "Another Name", "C2H6O", 46])
    sheet.append([4, "Formula Missing", None, 100])
    workbook.save(path)


def test_enfsi_entry_list_aggregate_crosscheck(tmp_path: Path) -> None:
    entry_list = tmp_path / "enfsi.xlsx"
    _write_enfsi(entry_list)
    dataset = tmp_path / "normalized.jsonl"
    rows = [
        {
            "compound_id": "cas:2759-28-6",
            "compound_name_raw": "1-Benzylpiperazine",
            "compound_name_companion": "1-Benzylpiperazine",
            "compound_name_canonical": "1-Benzylpiperazine",
            "molecular_formula": "C11H16N2",
            "molecular_weight": 176.0,
        },
        {
            "compound_id": "name:another",
            "compound_name_raw": "Another Name",
            "compound_name_companion": None,
            "compound_name_canonical": "Another Name",
            "molecular_formula": "C2H6O",
            "molecular_weight": 46.0,
        },
    ]
    dataset.write_text("\n".join(json.dumps(row) for row in rows) + "\n", encoding="utf-8")

    entries = read_enfsi_entry_list(entry_list)
    report = audit_enfsi_entry_list(entry_list, dataset)

    assert len(entries) == 4
    assert report.enfsi_entries == 4
    assert report.enfsi_unique_names == 3
    assert report.enfsi_missing_formula == 1
    assert report.enfsi_exact_duplicate_rows == 1
    assert report.high_confidence_matched_entries == 3
    assert report.high_confidence_matched_unique_swg_compounds == 2
    assert report.unmatched_entries == 1


def _mini_jcamp(title: str, name: str, peaks: list[tuple[int, int]]) -> str:
    peak_lines = "\n".join(f"{mz}, {intensity}" for mz, intensity in peaks)
    return (
        f"##TITLE={title}\n"
        "##JCAMP-DX=4.24\n"
        "##DATA TYPE=MASS SPECTRUM\n"
        f"##CAS NAME={name}\n"
        f"##NPOINTS={len(peaks)}\n"
        "##XYDATA=(XY..XY)\n"
        f"{peak_lines}\n"
        "##END=\n"
    )


def test_enfsi_spectral_overlap_detects_rescaled_near_duplicate() -> None:
    swgdrug = parse_jcamp_text(
        _mini_jcamp("Library Entry 1", "Compound A", [(50, 1000), (75, 500), (100, 100)]),
        source="swgdrug_eims",
        source_version="3.14",
    )
    enfsi = parse_jcamp_text(
        _mini_jcamp("Library Entry 1", "Compound A", [(50, 100), (75, 50), (100, 10)])
        + _mini_jcamp("Library Entry 2", "Compound B", [(40, 100), (60, 100)]),
        source="enfsi_dwg_eims",
        source_version="2025.11",
    )

    report, details = audit_enfsi_spectral_overlap(enfsi, swgdrug)

    assert report.enfsi_spectra == 2
    assert report.enfsi_duplicate_spectrum_records == 0
    assert report.exact_peak_list_match_records == 0
    assert report.identical_mz_support_match_records == 1
    assert report.near_identical_best_match_records == 1
    assert report.spectrum_disjoint_candidate_records == 1
    assert details[0].near_identical is True
    assert details[0].identical_mz_support is True
    assert details[0].cosine_similarity == 1.0
    assert details[1].near_identical is False


def test_read_jcamp_source_reads_single_member_zip_without_extracting(tmp_path: Path) -> None:
    archive_path = tmp_path / "enfsi.zip"
    payload = _mini_jcamp("Library Entry 1", "Compound A", [(50, 100), (75, 50)])
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.writestr("2025ENFSI.JDX", payload)

    records, provenance = read_jcamp_source(
        archive_path, source="enfsi_dwg_eims", source_version="2025.11"
    )

    assert len(records) == 1
    assert provenance["member"] == "2025ENFSI.JDX"
    assert len(str(provenance["archive_sha256"])) == 64
    assert len(str(provenance["member_sha256"])) == 64
    assert not (tmp_path / "2025ENFSI.JDX").exists()


def test_private_overlap_manifest_contains_ids_scores_but_not_names(tmp_path: Path) -> None:
    swgdrug = parse_jcamp_text(
        _mini_jcamp("Library Entry 1", "Secret SWG name", [(50, 1000), (75, 500)]),
        source="swgdrug_eims",
    )
    enfsi = parse_jcamp_text(
        _mini_jcamp("Library Entry 9", "Secret ENFSI name", [(50, 100), (75, 50)]),
        source="enfsi_dwg_eims",
    )
    _, details = audit_enfsi_spectral_overlap(enfsi, swgdrug)
    output = tmp_path / "private.jsonl"
    write_private_match_manifest(output, details)
    text = output.read_text(encoding="utf-8")

    assert "Library Entry 9" in text
    assert "Secret ENFSI name" not in text
    assert "Secret SWG name" not in text
