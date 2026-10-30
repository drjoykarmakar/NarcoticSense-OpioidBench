"""Conservative parsers for EI-MS source formats."""

from __future__ import annotations

import hashlib
import re
from collections.abc import Iterable
from datetime import date
from pathlib import Path

from opioidbench.data.schema import RawSpectrum, SourceProvenance, SpectrumPeak

_METADATA_RE = re.compile(r"^##\s*([^=]+?)\s*=\s*(.*)$")
_PAIR_RE = re.compile(
    r"(?P<mz>[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)\s*[, ]\s*"
    r"(?P<intensity>[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)"
)


class SpectrumParseError(ValueError):
    """Raised when a source record cannot be parsed safely."""


def _canonical_key(key: str) -> str:
    return re.sub(r"[^A-Z0-9]+", "_", key.strip().upper()).strip("_")


def _parse_peak_line(line: str) -> list[SpectrumPeak]:
    peaks: list[SpectrumPeak] = []
    for match in _PAIR_RE.finditer(line):
        peaks.append(SpectrumPeak(float(match.group("mz")), float(match.group("intensity"))))
    return peaks


def _record_id(metadata: dict[str, str], ordinal: int) -> str:
    # Prefer an explicit library record identifier.  In SWGDRUG HPJ, TITLE is
    # the stable "Library Entry N" key while CAS identifies the compound and
    # can repeat across replicate/reference spectra.
    for key in ("NIST_NO", "TITLE", "CAS_REGISTRY_NO", "CAS", "NAME"):
        value = metadata.get(key)
        if value:
            return value
    return f"record-{ordinal:05d}"


def _compound_name(metadata: dict[str, str]) -> str | None:
    """Return the source-provided chemical name without synonym inference.

    SWGDRUG HPJ files use ``TITLE`` for the library entry identifier and
    ``CAS NAME`` for the chemical name. Generic JCAMP files more commonly use
    ``TITLE``/``NAME`` directly, so those remain conservative fallbacks.
    """

    for key in ("CAS_NAME", "NAME", "TITLE"):
        value = metadata.get(key)
        if value and value.strip():
            return value.strip()
    return None


def _spectrum_id(source: str, record_id: str, peaks: Iterable[SpectrumPeak]) -> str:
    payload = [source, record_id]
    payload.extend(f"{peak.mz:.8g}:{peak.intensity:.8g}" for peak in peaks)
    digest = hashlib.sha256("|".join(payload).encode("utf-8")).hexdigest()[:20]
    return f"{source}:{digest}"


def parse_jcamp_text(
    text: str,
    *,
    source: str,
    source_version: str | None = None,
    retrieval_date: date | None = None,
    source_url: str | None = None,
    license_name: str | None = None,
    source_file: str | None = None,
    source_file_sha256: str | None = None,
) -> list[RawSpectrum]:
    """Parse one or more JCAMP-style EI-MS records.

    The parser intentionally preserves peak order and values. It supports common
    `##PEAK TABLE`/`##PEAKTABLE` and `##XYDATA` sections with explicit x,y pairs.
    Packed/differential JCAMP encodings are rejected rather than guessed.
    """

    records: list[RawSpectrum] = []
    metadata: dict[str, str] = {}
    peaks: list[SpectrumPeak] = []
    in_peak_section = False
    saw_record_content = False

    def finalize() -> None:
        nonlocal metadata, peaks, in_peak_section, saw_record_content
        if not saw_record_content:
            return
        ordinal = len(records) + 1
        rid = _record_id(metadata, ordinal)
        records.append(
            RawSpectrum(
                spectrum_id=_spectrum_id(source, rid, peaks),
                compound_name_raw=_compound_name(metadata),
                peaks=tuple(peaks),
                provenance=SourceProvenance(
                    source=source,
                    source_version=source_version,
                    source_record_id=rid,
                    retrieval_date=retrieval_date,
                    source_url=source_url,
                    license=license_name,
                    source_file=source_file,
                    source_file_sha256=source_file_sha256,
                ),
                metadata=dict(metadata),
            )
        )
        metadata = {}
        peaks = []
        in_peak_section = False
        saw_record_content = False

    for line_number, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.upper().startswith("##END"):
            finalize()
            continue
        meta_match = _METADATA_RE.match(line)
        if meta_match:
            key = _canonical_key(meta_match.group(1))
            value = meta_match.group(2).strip()
            # SWGDRUG HPJ is a concatenated JCAMP stream without ##END=.  A
            # repeated TITLE therefore marks the start of the next record.
            if key == "TITLE" and saw_record_content and "TITLE" in metadata:
                finalize()
            saw_record_content = True
            metadata[key] = value
            in_peak_section = key in {"PEAK_TABLE", "PEAKTABLE", "XYDATA", "XYPOINTS"}
            continue
        if line.startswith("##"):
            saw_record_content = True
            continue
        if in_peak_section:
            parsed = _parse_peak_line(line)
            if not parsed:
                raise SpectrumParseError(
                    f"Unsupported or malformed JCAMP peak data at line {line_number}: {line!r}"
                )
            peaks.extend(parsed)
            saw_record_content = True

    finalize()
    if not records:
        raise SpectrumParseError("No JCAMP records found")
    return records


def parse_jcamp_file(
    path: Path,
    *,
    source: str,
    source_version: str | None = None,
    retrieval_date: date | None = None,
    source_url: str | None = None,
    license_name: str | None = None,
    source_file_sha256: str | None = None,
) -> list[RawSpectrum]:
    """Parse a JCAMP file as UTF-8 with conservative replacement disabled."""
    text = path.read_text(encoding="utf-8")
    return parse_jcamp_text(
        text,
        source=source,
        source_version=source_version,
        retrieval_date=retrieval_date,
        source_url=source_url,
        license_name=license_name,
        source_file=str(path),
        source_file_sha256=source_file_sha256,
    )
