from datetime import date
from pathlib import Path

from opioidbench.data.parsers import parse_jcamp_file
from opioidbench.data.validation import normalize_spectrum, validate_dataset

FIXTURES = Path(__file__).parent / "fixtures"


def _parse(name: str):
    return parse_jcamp_file(
        FIXTURES / name,
        source="synthetic_fixture",
        source_version="1",
        retrieval_date=date(2026, 9, 29),
    )


def test_valid_fixture_passes_qc_and_preserves_missing_metadata() -> None:
    validated, report = validate_dataset(_parse("synthetic_eims.jdx"))
    assert report.number_spectra == 2
    assert report.number_invalid == 0
    normalized = normalize_spectrum(validated[0])
    assert normalized.compound_id.startswith("name:")
    assert normalized.cas is None
    assert normalized.canonical_smiles is None
    assert normalized.opioid_label is None
    assert len(normalized.peaks) == 4


def test_swgdrug_style_placeholder_cas_is_preserved_as_qc_but_normalized_missing() -> None:
    validated, report = validate_dataset(_parse("synthetic_swgdrug_layout.hpj"))
    assert report.issue_counts["placeholder_cas"] == 2
    normalized = normalize_spectrum(validated[0])
    assert normalized.cas is None
    assert normalized.molecular_weight == 46.0


def test_invalid_fixture_reports_core_peak_problems() -> None:
    _, report = validate_dataset(_parse("invalid_eims.jdx"))
    assert report.number_invalid == 1
    assert report.issue_counts["missing_compound_name"] == 1
    assert report.issue_counts["invalid_mz"] == 1
    assert report.issue_counts["negative_intensity"] == 1
    assert report.issue_counts["duplicate_peaks"] == 1


def test_duplicate_spectra_are_flagged() -> None:
    records = _parse("synthetic_eims.jdx")
    _, report = validate_dataset([records[0], records[0]])
    assert report.issue_counts["duplicate_spectrum"] == 1
    assert report.issue_counts["duplicate_spectrum_id"] == 1
