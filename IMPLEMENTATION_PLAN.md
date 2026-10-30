# First Milestone Implementation Plan

## Scope

Implement only the publication-grade foundation needed to demonstrate:

`raw EI-MS spectrum -> parse -> schema validation -> QC -> normalized record`

No cheminformatics enrichment, taxonomy curation, benchmark splitting, ML, OOD, uncertainty, or CNN code is included in this milestone.

## Files to create

- Project/reproducibility: `pyproject.toml`, `environment.yml`, `Makefile`, `.gitignore`, `.pre-commit-config.yaml`, CI workflow.
- Governance/docs: `README.md`, `LICENSE`, `CITATION.cff`, contributing/security/code-of-conduct/changelog, docs skeleton, data README, dataset-card skeleton.
- Source layer: `metadata/source_registry.yaml`, source-adapter base class, SWGDRUG adapter, safe downloader/importer with checksums.
- EI-MS data layer: dataclass schemas, JCAMP parser, validation/QC, normalized JSONL output.
- Interface: CLI and thin scripts for source inspection/import, validation, and synthetic demo.
- Tests: synthetic fixtures covering parsing, schema invariants, QC, source-policy enforcement, checksum verification, and CLI demo.

## Dependencies

Runtime is intentionally small: Python >=3.11 and PyYAML. Development tooling uses pytest and ruff. No RDKit, scikit-learn, or PyTorch is required for the first milestone.

## Data-source strategy

All source-specific behavior is behind a `DataSource` interface. The SWGDRUG registry entry records version/URL/retrieval metadata but does not assume a redistribution license. Automated download is blocked unless the registry explicitly marks it allowed. Locally and legally obtained source files can be imported with SHA-256 provenance.

## Test strategy

Tests use only tiny synthetic JCAMP files clearly marked synthetic. They verify parser behavior, deterministic IDs/hashes, duplicate/invalid peak detection, dataset-level duplicate detection, source-registry loading, download-policy enforcement, checksum mismatch handling, and an end-to-end CLI demo.

## Licensing/data-access concerns

The SWGDRUG website makes the library downloadable and provides a disclaimer, but an explicit redistribution license was not identified during this milestone. Raw SWGDRUG data therefore must not be committed by default. The registry marks redistribution as `unresolved`, and the code requires local import unless a future documented review changes the policy.
