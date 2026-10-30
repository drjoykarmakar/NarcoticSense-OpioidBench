from datetime import date
from pathlib import Path

import pytest

from opioidbench.data.parsers import SpectrumParseError, parse_jcamp_file, parse_jcamp_text

FIXTURES = Path(__file__).parent / "fixtures"


def test_parse_synthetic_jcamp_records() -> None:
    records = parse_jcamp_file(
        FIXTURES / "synthetic_eims.jdx",
        source="synthetic_fixture",
        source_version="1",
        retrieval_date=date(2026, 9, 29),
    )
    assert len(records) == 2
    assert records[0].compound_name_raw == "Synthetic Example A"
    assert records[0].peaks[1].mz == 57
    assert records[0].peaks[1].intensity == 100
    assert records[0].metadata["ORIGIN"].startswith("OpioidBench synthetic fixture")
    assert records[0].spectrum_id.startswith("synthetic_fixture:")


def test_parser_rejects_malformed_peak_encoding() -> None:
    with pytest.raises(SpectrumParseError, match="malformed"):
        parse_jcamp_text(
            "##TITLE=X\n##PEAK TABLE=(XY..XY)\nthis-is-not-a-pair\n##END=",
            source="test",
        )


def test_parse_swgdrug_style_concatenated_records_without_end_markers() -> None:
    records = parse_jcamp_file(
        FIXTURES / "synthetic_swgdrug_layout.hpj",
        source="synthetic_fixture",
        source_version="1",
    )
    assert len(records) == 2
    assert records[0].provenance.source_record_id == "Library Entry 1001"
    assert records[0].compound_name_raw == "Synthetic Compound Alpha"
    assert len(records[0].peaks) == 3
    assert records[1].provenance.source_record_id == "Library Entry 1002"
    assert records[1].compound_name_raw == "Synthetic Compound Beta"
