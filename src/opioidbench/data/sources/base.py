"""Source registry and adapter abstraction."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


@dataclass(frozen=True, slots=True)
class SourceSpec:
    key: str
    display_name: str
    modality: str
    source_url: str | None
    version: str | None
    version_date: str | None
    retrieval_date: str | None
    file_format: str | None
    license: str | None
    terms_url: str | None
    redistribution_status: str
    automatic_download_allowed: bool
    citation_requirements: str | None
    notes: str | None
    observed_archive_sha256: str | None = None
    observed_member_sha256: str | None = None
    observed_companion_sha256: str | None = None
    observed_exact_mass_archive_sha256: str | None = None
    observed_exact_mass_workbook_sha256: str | None = None
    observed_exact_mass_pcdl_sha256: str | None = None
    observed_agilent_archive_sha256: str | None = None
    observed_agilent_header_sha256: str | None = None
    observed_nist_archive_sha256: str | None = None
    observed_nist_structure_sha256: str | None = None
    observed_shimadzu_archive_sha256: str | None = None
    observed_changelog_sha256: str | None = None
    observed_entry_list_sha256: str | None = None


class SourceRegistry:
    """Load auditable source acquisition policy from YAML."""

    def __init__(self, path: Path):
        payload = yaml.safe_load(path.read_text(encoding="utf-8"))
        if not isinstance(payload, dict) or not isinstance(payload.get("sources"), dict):
            raise ValueError("Source registry must contain a 'sources' mapping")
        self.schema_version = int(payload.get("schema_version", 1))
        self._sources: dict[str, SourceSpec] = {}
        for key, raw in payload["sources"].items():
            if not isinstance(raw, dict):
                raise ValueError(f"Source entry {key!r} must be a mapping")
            self._sources[key] = SourceSpec(key=key, **raw)

    def get(self, key: str) -> SourceSpec:
        try:
            return self._sources[key]
        except KeyError as exc:
            raise KeyError(f"Unknown source {key!r}") from exc

    def list(self) -> tuple[SourceSpec, ...]:
        return tuple(self._sources[key] for key in sorted(self._sources))


class DataSource(ABC):
    """Minimal adapter interface separating acquisition from ML code."""

    def __init__(self, spec: SourceSpec):
        self.spec = spec

    @abstractmethod
    def supported_suffixes(self) -> tuple[str, ...]:
        """Return source file suffixes accepted by this adapter."""

    def describe(self) -> dict[str, Any]:
        return {
            "key": self.spec.key,
            "display_name": self.spec.display_name,
            "modality": self.spec.modality,
            "version": self.spec.version,
            "source_url": self.spec.source_url,
            "license": self.spec.license,
            "redistribution_status": self.spec.redistribution_status,
            "automatic_download_allowed": self.spec.automatic_download_allowed,
        }
