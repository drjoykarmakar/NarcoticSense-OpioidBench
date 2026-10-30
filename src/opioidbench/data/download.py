"""Safe acquisition/import helpers with explicit policy and checksums."""

from __future__ import annotations

import hashlib
import shutil
import urllib.request
import zipfile
from pathlib import Path

from opioidbench.data.sources.base import SourceSpec


class DataAccessPolicyError(RuntimeError):
    """Raised when source policy disallows an automated acquisition action."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_sha256(path: Path, expected: str) -> str:
    actual = sha256_file(path)
    if actual.lower() != expected.lower():
        raise ValueError(f"SHA-256 mismatch for {path}: expected {expected}, got {actual}")
    return actual


def import_local_file(
    source_path: Path,
    destination_dir: Path,
    *,
    expected_sha256: str | None = None,
) -> tuple[Path, str]:
    """Copy a legally obtained local file into the ignored raw-data area."""
    if not source_path.is_file():
        raise FileNotFoundError(source_path)
    actual = sha256_file(source_path)
    if expected_sha256 and actual.lower() != expected_sha256.lower():
        raise ValueError(
            f"SHA-256 mismatch for {source_path}: expected {expected_sha256}, got {actual}"
        )
    destination_dir.mkdir(parents=True, exist_ok=True)
    destination = destination_dir / source_path.name
    shutil.copy2(source_path, destination)
    return destination, actual


def import_local_zip(
    archive_path: Path,
    destination_dir: Path,
    *,
    member_name: str | None = None,
    expected_archive_sha256: str | None = None,
) -> tuple[Path, str, str, str]:
    """Import one spectrum-library member from a locally supplied ZIP archive."""
    if not archive_path.is_file():
        raise FileNotFoundError(archive_path)
    archive_sha256 = sha256_file(archive_path)
    if expected_archive_sha256 and archive_sha256.lower() != expected_archive_sha256.lower():
        raise ValueError(
            f"SHA-256 mismatch for {archive_path}: expected {expected_archive_sha256}, "
            f"got {archive_sha256}"
        )
    allowed_suffixes = {".hpj", ".jdx", ".dx", ".txt"}
    with zipfile.ZipFile(archive_path) as archive:
        candidates = [
            info
            for info in archive.infolist()
            if not info.is_dir() and Path(info.filename).suffix.lower() in allowed_suffixes
        ]
        if member_name:
            candidates = [info for info in candidates if info.filename == member_name]
        if len(candidates) != 1:
            names = [info.filename for info in candidates]
            raise ValueError(
                "ZIP import requires exactly one supported spectrum member; "
                f"found {len(candidates)}: {names}"
            )
        info = candidates[0]
        destination_dir.mkdir(parents=True, exist_ok=True)
        destination = destination_dir / Path(info.filename).name
        with archive.open(info) as source, destination.open("wb") as target:
            shutil.copyfileobj(source, target)
    member_sha256 = sha256_file(destination)
    return destination, archive_sha256, member_sha256, info.filename


def download_source(spec: SourceSpec, destination: Path) -> tuple[Path, str]:
    """Download only sources whose registry policy explicitly permits it."""
    if not spec.automatic_download_allowed:
        raise DataAccessPolicyError(
            f"Automatic download is disabled for {spec.key}. Obtain the source legally, then "
            f"run: opioidbench data import --source {spec.key} --input /path/to/file "
            f"--destination data/raw/{spec.key}. Review metadata/source_registry.yaml for "
            "the documented access and redistribution policy."
        )
    if not spec.source_url:
        raise DataAccessPolicyError(f"Source {spec.key} has no downloadable URL")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(spec.source_url, timeout=60) as response:
        destination.write_bytes(response.read())
    return destination, sha256_file(destination)
