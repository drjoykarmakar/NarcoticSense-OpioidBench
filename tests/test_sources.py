from pathlib import Path

import pytest

from opioidbench.data.download import DataAccessPolicyError, download_source, import_local_file
from opioidbench.data.sources.base import SourceRegistry

ROOT = Path(__file__).parents[1]


def test_registry_loads_swgdrug_with_conservative_policy() -> None:
    registry = SourceRegistry(ROOT / "metadata" / "source_registry.yaml")
    spec = registry.get("swgdrug_eims")
    assert spec.version == "3.14"
    assert spec.redistribution_status == "unresolved"
    assert spec.automatic_download_allowed is False


def test_swgdrug_automatic_download_is_blocked(tmp_path: Path) -> None:
    spec = SourceRegistry(ROOT / "metadata" / "source_registry.yaml").get("swgdrug_eims")
    with pytest.raises(DataAccessPolicyError, match="Automatic download is disabled"):
        download_source(spec, tmp_path / "swgdrug.dat")


def test_local_import_records_checksum_and_verifies_expected(tmp_path: Path) -> None:
    source = tmp_path / "source.jdx"
    source.write_text("synthetic", encoding="utf-8")
    imported, digest = import_local_file(source, tmp_path / "raw")
    assert imported.read_text(encoding="utf-8") == "synthetic"
    assert len(digest) == 64
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        import_local_file(source, tmp_path / "other", expected_sha256="0" * 64)


def test_packaged_registry_matches_reviewable_registry() -> None:
    assert (ROOT / "metadata" / "source_registry.yaml").read_bytes() == (
        ROOT / "src" / "opioidbench" / "resources" / "source_registry.yaml"
    ).read_bytes()


def test_registry_marks_enfsi_as_restricted_local_only() -> None:
    registry = SourceRegistry(ROOT / "metadata" / "source_registry.yaml")
    spec = registry.get("enfsi_dwg_eims")
    assert spec.version == "2025.11"
    assert spec.redistribution_status == "restricted"
    assert spec.automatic_download_allowed is False
    assert spec.observed_entry_list_sha256 == (
        "55cafbf460edc96e88c5a499a95d04c71b1ba34c5d144cdc8af1acc26bee471b"
    )
    assert spec.observed_archive_sha256 == (
        "90f5305c9bb431c2e562dd283a05e3fe001be57168d15a1a546f803f55de8946"
    )
    assert spec.observed_member_sha256 == (
        "f7f96fc24dda46694e9383aea2b1d98025a6cd03f584479503a41213ea450146"
    )


def test_registry_tracks_chebi_bulk_snapshot() -> None:
    registry = SourceRegistry(ROOT / "metadata" / "source_registry.yaml")
    spec = registry.get("chebi_full")
    assert spec.modality == "chemical_metadata"
    assert spec.license == "CC BY 4.0"
    assert spec.redistribution_status == "allowed_with_attribution"
    assert spec.observed_archive_sha256 == (
        "9e180253b1d6b57c197a06e5ea97c500488935c85a2df9ffcea6339e56ab3451"
    )
