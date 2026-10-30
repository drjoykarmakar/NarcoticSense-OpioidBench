"""Command-line interface for the first OpioidBench milestone."""

from __future__ import annotations

import argparse
import json
import logging
from datetime import date
from importlib.resources import files
from pathlib import Path

from opioidbench.chemistry.enrichment import enrich_compound_table
from opioidbench.chemistry.evidence import curate_taxonomy_from_evidence
from opioidbench.chemistry.resolution import resolve_compound_sources
from opioidbench.chemistry.review import prioritize_taxonomy_review
from opioidbench.data.build import build_dataset
from opioidbench.data.companion import (
    align_swgdrug_agilent_header,
    align_swgdrug_compound_list,
    align_swgdrug_exact_mass_list,
    read_swgdrug_agilent_header_archive,
    read_swgdrug_compound_list,
    read_swgdrug_exact_mass_list,
)
from opioidbench.data.download import (
    DataAccessPolicyError,
    download_source,
    import_local_file,
    import_local_zip,
    sha256_file,
)
from opioidbench.data.enfsi import audit_enfsi_entry_list
from opioidbench.data.external import (
    audit_enfsi_external_candidates,
    write_private_external_candidate_manifest,
)
from opioidbench.data.nist import (
    align_swgdrug_nist_structures,
    read_swgdrug_nist_structure_archive,
)
from opioidbench.data.overlap import (
    audit_enfsi_spectral_overlap,
    read_jcamp_source,
    write_private_match_manifest,
)
from opioidbench.data.parsers import SpectrumParseError, parse_jcamp_file
from opioidbench.data.sources.base import SourceRegistry
from opioidbench.data.validation import normalize_spectrum, validate_dataset
from opioidbench.experiments.reproducibility import write_json

LOGGER = logging.getLogger("opioidbench")
PROJECT_REGISTRY = Path(__file__).resolve().parents[2] / "metadata" / "source_registry.yaml"
PROJECT_EVIDENCE = Path(__file__).resolve().parents[2] / "metadata" / "opioid_evidence.csv"
PROJECT_EXTERNAL_STRUCTURES = (
    Path(__file__).resolve().parents[2] / "metadata" / "external_structure_evidence.csv"
)
PROJECT_TAXONOMY_EXCLUSIONS = (
    Path(__file__).resolve().parents[2] / "metadata" / "taxonomy_exclusions.csv"
)


def _default_registry_path() -> Path:
    if PROJECT_REGISTRY.is_file():
        return PROJECT_REGISTRY
    resource = files("opioidbench").joinpath("resources/source_registry.yaml")
    return Path(str(resource))


def _default_evidence_path() -> Path:
    if PROJECT_EVIDENCE.is_file():
        return PROJECT_EVIDENCE
    resource = files("opioidbench").joinpath("resources/opioid_evidence.csv")
    return Path(str(resource))


def _default_external_structures_path() -> Path:
    if PROJECT_EXTERNAL_STRUCTURES.is_file():
        return PROJECT_EXTERNAL_STRUCTURES
    resource = files("opioidbench").joinpath("resources/external_structure_evidence.csv")
    return Path(str(resource))


def _default_taxonomy_exclusions_path() -> Path:
    if PROJECT_TAXONOMY_EXCLUSIONS.is_file():
        return PROJECT_TAXONOMY_EXCLUSIONS
    resource = files("opioidbench").joinpath("resources/taxonomy_exclusions.csv")
    return Path(str(resource))


def _configure_logging(verbose: bool = False) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(levelname)s %(name)s %(message)s",
    )


def _registry(path: str | None) -> SourceRegistry:
    return SourceRegistry(Path(path) if path else _default_registry_path())


def _parse_with_source(input_path: Path, source_key: str, registry_path: str | None):
    registry = _registry(registry_path)
    spec = registry.get(source_key)
    digest = sha256_file(input_path)
    sidecar = input_path.with_suffix(input_path.suffix + ".provenance.json")
    retrieval = None
    if sidecar.is_file():
        provenance = json.loads(sidecar.read_text(encoding="utf-8"))
        if provenance.get("retrieval_date"):
            retrieval = date.fromisoformat(provenance["retrieval_date"])
    elif spec.retrieval_date:
        retrieval = date.fromisoformat(spec.retrieval_date)
    records = parse_jcamp_file(
        input_path,
        source=spec.key,
        source_version=spec.version,
        retrieval_date=retrieval,
        source_url=spec.source_url,
        license_name=spec.license,
        source_file_sha256=digest,
    )
    return spec, records, digest


def cmd_sources(args: argparse.Namespace) -> int:
    for spec in _registry(args.registry).list():
        print(
            json.dumps(
                {
                    "key": spec.key,
                    "display_name": spec.display_name,
                    "version": spec.version,
                    "redistribution_status": spec.redistribution_status,
                    "automatic_download_allowed": spec.automatic_download_allowed,
                },
                sort_keys=True,
            )
        )
    return 0


def cmd_download(args: argparse.Namespace) -> int:
    spec = _registry(args.registry).get(args.source)
    destination, digest = download_source(spec, Path(args.destination))
    LOGGER.info("Downloaded %s with SHA-256 %s", destination, digest)
    return 0


def cmd_import(args: argparse.Namespace) -> int:
    spec = _registry(args.registry).get(args.source)
    input_path = Path(args.input)
    archive_digest = None
    archive_member = None
    if input_path.suffix.lower() == ".zip":
        destination, archive_digest, digest, archive_member = import_local_zip(
            input_path,
            Path(args.destination),
            member_name=args.member,
            expected_archive_sha256=args.sha256,
        )
    else:
        destination, digest = import_local_file(
            input_path, Path(args.destination), expected_sha256=args.sha256
        )
    manifest = {
        "source": spec.key,
        "source_version": spec.version,
        "source_url": spec.source_url,
        "license": spec.license,
        "redistribution_status": spec.redistribution_status,
        "local_file": str(destination),
        "sha256": digest,
        "source_archive": input_path.name if archive_digest else None,
        "source_archive_sha256": archive_digest,
        "source_archive_member": archive_member,
        "import_date": date.today().isoformat(),
        "retrieval_date": args.retrieval_date,
    }
    manifest_path = destination.with_suffix(destination.suffix + ".provenance.json")
    write_json(manifest_path, manifest)
    LOGGER.info("Imported %s with SHA-256 %s", destination, digest)
    return 0


def cmd_validate(args: argparse.Namespace) -> int:
    _, records, digest = _parse_with_source(Path(args.input), args.source, args.registry)
    _, report = validate_dataset(records)
    payload = report.to_dict()
    payload["source_file_sha256"] = digest
    if args.output:
        write_json(Path(args.output), payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 1 if report.number_invalid else 0


def cmd_build(args: argparse.Namespace) -> int:
    input_path = Path(args.input)
    _, records, digest = _parse_with_source(input_path, args.source, args.registry)
    companion_path = Path(args.compound_list) if args.compound_list else None
    companion_digest = None
    companion_alignment = None
    if companion_path:
        entries, _ = read_swgdrug_compound_list(companion_path)
        records, alignment = align_swgdrug_compound_list(records, entries)
        companion_digest = sha256_file(companion_path)
        companion_alignment = alignment.to_dict()

    exact_mass_path = Path(args.exact_mass_list) if args.exact_mass_list else None
    exact_mass_digest = None
    exact_mass_alignment = None
    if exact_mass_path:
        exact_entries, exact_version = read_swgdrug_exact_mass_list(exact_mass_path)
        records, exact_alignment = align_swgdrug_exact_mass_list(records, exact_entries)
        exact_mass_digest = sha256_file(exact_mass_path)
        exact_mass_alignment = exact_alignment.to_dict()
        exact_mass_alignment["source_version_marker"] = exact_version

    agilent_archive = Path(args.agilent_archive) if args.agilent_archive else None
    agilent_archive_digest = None
    agilent_header_alignment = None
    if agilent_archive:
        agilent_entries, header_member, header_digest = read_swgdrug_agilent_header_archive(
            agilent_archive
        )
        records, agilent_alignment = align_swgdrug_agilent_header(
            records,
            agilent_entries,
            archive_member=header_member,
            header_sha256=header_digest,
        )
        agilent_archive_digest = sha256_file(agilent_archive)
        agilent_header_alignment = agilent_alignment.to_dict()

    nist_archive = Path(args.nist_archive) if args.nist_archive else None
    nist_archive_digest = None
    nist_structure_member = None
    nist_structure_member_digest = None
    nist_structure_alignment = None
    nist_structure_review = None
    if nist_archive:
        nist_entries, nist_structure_member, nist_structure_member_digest = (
            read_swgdrug_nist_structure_archive(nist_archive)
        )
        records, nist_alignment, nist_structure_review = align_swgdrug_nist_structures(
            records,
            nist_entries,
            archive_member=nist_structure_member,
            member_sha256=nist_structure_member_digest,
        )
        nist_archive_digest = sha256_file(nist_archive)
        nist_structure_alignment = nist_alignment.to_dict()

    archive_path = Path(args.source_archive) if args.source_archive else None
    archive_digest = sha256_file(archive_path) if archive_path else None
    result = build_dataset(
        records,
        output_dir=Path(args.output),
        dataset_version=args.dataset_version,
        source_file=input_path,
        source_file_sha256=digest,
        taxonomy_path=Path(args.taxonomy) if args.taxonomy else None,
        source_archive=archive_path,
        source_archive_sha256=archive_digest,
        companion_file=companion_path,
        companion_file_sha256=companion_digest,
        companion_alignment=companion_alignment,
        exact_mass_file=exact_mass_path,
        exact_mass_file_sha256=exact_mass_digest,
        exact_mass_alignment=exact_mass_alignment,
        agilent_archive=agilent_archive,
        agilent_archive_sha256=agilent_archive_digest,
        agilent_header_alignment=agilent_header_alignment,
        nist_archive=nist_archive,
        nist_archive_sha256=nist_archive_digest,
        nist_structure_member=nist_structure_member,
        nist_structure_member_sha256=nist_structure_member_digest,
        nist_structure_alignment=nist_structure_alignment,
        nist_structure_review=nist_structure_review,
    )
    print(json.dumps({"audit": result["audit"], "manifest": result["manifest"]}, indent=2))
    return 1 if result["audit"]["number_invalid_spectra"] else 0


def cmd_demo(args: argparse.Namespace) -> int:
    input_path = Path(args.fixture)
    output_dir = Path(args.output)
    spec, records, digest = _parse_with_source(input_path, "synthetic_fixture", args.registry)
    validated, report = validate_dataset(records)
    output_dir.mkdir(parents=True, exist_ok=True)
    write_json(output_dir / "data_qc.json", {**report.to_dict(), "source_file_sha256": digest})

    normalized_path = output_dir / "normalized_records.jsonl"
    with normalized_path.open("w", encoding="utf-8") as handle:
        for item in validated:
            handle.write(json.dumps(normalize_spectrum(item).to_dict(), sort_keys=True) + "\n")

    manifest = {
        "stage": "first_milestone_demo",
        "source": spec.key,
        "source_version": spec.version,
        "source_file": str(input_path),
        "source_file_sha256": digest,
        "number_spectra": len(records),
        "number_valid": report.number_valid,
        "number_invalid": report.number_invalid,
        "outputs": [str(output_dir / "data_qc.json"), str(normalized_path)],
        "scientific_result": False,
    }
    write_json(output_dir / "manifest.json", manifest)
    LOGGER.info(
        "Parsed %d spectrum(s); valid=%d invalid=%d; output=%s",
        len(records),
        report.number_valid,
        report.number_invalid,
        output_dir,
    )
    return 1 if report.number_invalid else 0


def cmd_data_audit_enfsi_overlap(args: argparse.Namespace) -> int:
    registry = _registry(args.registry)
    enfsi_spec = registry.get("enfsi_dwg_eims")
    swgdrug_spec = registry.get("swgdrug_eims")
    enfsi_path = Path(args.enfsi)
    swgdrug_path = Path(args.swgdrug)
    enfsi_records, enfsi_provenance = read_jcamp_source(
        enfsi_path, source=enfsi_spec.key, source_version=enfsi_spec.version or "unknown"
    )
    swgdrug_records, swgdrug_provenance = read_jcamp_source(
        swgdrug_path, source=swgdrug_spec.key, source_version=swgdrug_spec.version or "unknown"
    )
    report, details = audit_enfsi_spectral_overlap(
        enfsi_records,
        swgdrug_records,
        near_identical_threshold=args.near_identical_threshold,
        entry_list_path=Path(args.entry_list) if args.entry_list else None,
        normalized_swgdrug_path=(Path(args.swgdrug_dataset) if args.swgdrug_dataset else None),
    )
    payload = report.to_dict()
    payload.update(
        {
            "enfsi_source_archive": (
                enfsi_path.name if enfsi_path.suffix.lower() == ".zip" else None
            ),
            "enfsi_source_archive_sha256": enfsi_provenance["archive_sha256"],
            "enfsi_source_member": enfsi_provenance["member"],
            "enfsi_source_member_sha256": enfsi_provenance["member_sha256"],
            "swgdrug_source_file": swgdrug_path.name,
            "swgdrug_source_sha256": swgdrug_provenance["member_sha256"],
            "entry_list_sha256": (sha256_file(Path(args.entry_list)) if args.entry_list else None),
            "restricted_source_data_in_aggregate_output": False,
            "interpretation": (
                "Near-identical is a conservative spectral leakage-risk flag, not proof of "
                "a shared physical acquisition. Spectrum-disjoint candidates are not "
                "necessarily compound-disjoint."
            ),
        }
    )
    if args.output:
        write_json(Path(args.output), payload)
    if args.private_output:
        write_private_match_manifest(Path(args.private_output), details)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


def cmd_data_audit_enfsi_external(args: argparse.Namespace) -> int:
    registry = _registry(args.registry)
    enfsi_spec = registry.get("enfsi_dwg_eims")
    swgdrug_spec = registry.get("swgdrug_eims")
    enfsi_path = Path(args.enfsi)
    swgdrug_path = Path(args.swgdrug)
    enfsi_records, enfsi_provenance = read_jcamp_source(
        enfsi_path, source=enfsi_spec.key, source_version=enfsi_spec.version or "unknown"
    )
    swgdrug_records, swgdrug_provenance = read_jcamp_source(
        swgdrug_path, source=swgdrug_spec.key, source_version=swgdrug_spec.version or "unknown"
    )
    report, private_rows = audit_enfsi_external_candidates(
        enfsi_records,
        swgdrug_records,
        entry_list_path=Path(args.entry_list),
        normalized_swgdrug_path=Path(args.swgdrug_dataset),
        resolved_swgdrug_compounds_path=Path(args.swgdrug_compounds),
        chebi_sdf_path=Path(args.chebi_sdf),
        supplemental_structures_path=(
            Path(args.supplemental_structures) if args.supplemental_structures else None
        ),
        taxonomy_evidence_path=(Path(args.taxonomy_evidence) if args.taxonomy_evidence else None),
        near_identical_threshold=args.near_identical_threshold,
    )
    report.update(
        {
            "enfsi_source_archive": (
                enfsi_path.name if enfsi_path.suffix.lower() == ".zip" else None
            ),
            "enfsi_source_archive_sha256": enfsi_provenance["archive_sha256"],
            "enfsi_source_member": enfsi_provenance["member"],
            "enfsi_source_member_sha256": enfsi_provenance["member_sha256"],
            "swgdrug_source_file": swgdrug_path.name,
            "swgdrug_source_sha256": swgdrug_provenance["member_sha256"],
            "entry_list_sha256": sha256_file(Path(args.entry_list)),
            "swgdrug_dataset_sha256": sha256_file(Path(args.swgdrug_dataset)),
            "swgdrug_compounds_sha256": sha256_file(Path(args.swgdrug_compounds)),
            "chebi_sdf_sha256": sha256_file(Path(args.chebi_sdf)),
            "supplemental_structures_sha256": (
                sha256_file(Path(args.supplemental_structures))
                if args.supplemental_structures
                else None
            ),
            "taxonomy_evidence_sha256": (
                sha256_file(Path(args.taxonomy_evidence)) if args.taxonomy_evidence else None
            ),
        }
    )
    if args.output:
        write_json(Path(args.output), report)
    if args.private_output:
        write_private_external_candidate_manifest(Path(args.private_output), private_rows)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


def cmd_metadata_enrich(args: argparse.Namespace) -> int:
    counts = enrich_compound_table(
        Path(args.input),
        output_dir=Path(args.output),
        cache_dir=Path(args.cache),
        allow_network=args.allow_network,
        max_records=args.max_records,
    )
    print(json.dumps(counts, indent=2, sort_keys=True))
    return 0 if counts["errors"] == 0 else 1


def cmd_metadata_resolve(args: argparse.Namespace) -> int:
    counts = resolve_compound_sources(
        Path(args.input),
        output_dir=Path(args.output),
        chebi_sdf=Path(args.chebi_sdf) if args.chebi_sdf else None,
        pubchem_cache_dir=Path(args.cache),
        allow_pubchem_network=args.allow_network,
        max_records=args.max_records,
    )
    print(json.dumps(counts, indent=2, sort_keys=True))
    return 0 if counts["errors"] == 0 else 1


def cmd_metadata_curate_taxonomy(args: argparse.Namespace) -> int:
    counts = curate_taxonomy_from_evidence(
        Path(args.input),
        Path(args.evidence) if args.evidence else _default_evidence_path(),
        taxonomy_output=Path(args.output),
        review_output=Path(args.review_output),
        exclusions_path=(
            Path(args.exclusions) if args.exclusions else _default_taxonomy_exclusions_path()
        ),
    )
    print(json.dumps(counts, indent=2, sort_keys=True))
    return 0 if counts["uncertain_conflict"] == 0 else 1


def cmd_metadata_prioritize_review(args: argparse.Namespace) -> int:
    counts = prioritize_taxonomy_review(
        Path(args.input),
        Path(args.taxonomy),
        output_path=Path(args.output),
        threshold=args.threshold,
    )
    print(json.dumps(counts, indent=2, sort_keys=True))
    return 0


def cmd_metadata_audit_enfsi(args: argparse.Namespace) -> int:
    report = audit_enfsi_entry_list(Path(args.entry_list), Path(args.dataset))
    payload = report.to_dict()
    payload["entry_list_sha256"] = sha256_file(Path(args.entry_list))
    payload["restricted_source_data_in_output"] = False
    if args.output:
        write_json(Path(args.output), payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


def cmd_benchmark_split(args: argparse.Namespace) -> int:
    from opioidbench.experiments.splits import generate_cv_splits

    report = generate_cv_splits(
        Path(args.dataset),
        Path(args.compounds),
        Path(args.taxonomy),
        output_dir=Path(args.output),
        protocol=args.protocol,
        n_splits=args.folds,
        seed=args.seed,
        near_identical_threshold=args.near_identical_threshold,
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


def cmd_benchmark_baselines(args: argparse.Namespace) -> int:
    from opioidbench.experiments.baselines import evaluate_baselines
    from opioidbench.experiments.splits import load_curated_samples, load_fold_assignments

    require_scaffold = args.protocol in {"scaffold", "analogue"}
    samples, _ = load_curated_samples(
        Path(args.dataset),
        Path(args.compounds),
        Path(args.taxonomy),
        require_scaffold=require_scaffold,
    )
    assignments = load_fold_assignments(Path(args.assignments))
    report = evaluate_baselines(
        samples,
        assignments,
        output_dir=Path(args.output),
        mz_max=args.mz_max,
        seed=args.seed,
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="opioidbench")
    parser.add_argument("--verbose", action="store_true")
    sub = parser.add_subparsers(dest="command", required=True)

    data = sub.add_parser("data", help="Data acquisition and validation")
    data_sub = data.add_subparsers(dest="data_command", required=True)

    sources = data_sub.add_parser("sources", help="List registered data sources")
    sources.add_argument("--registry")
    sources.set_defaults(func=cmd_sources)

    downloader = data_sub.add_parser("download", help="Download a source only when policy permits")
    downloader.add_argument("--source", required=True)
    downloader.add_argument("--destination", required=True)
    downloader.add_argument("--registry")
    downloader.set_defaults(func=cmd_download)

    importer = data_sub.add_parser("import", help="Import a legally obtained local source file")
    importer.add_argument("--source", required=True)
    importer.add_argument("--input", required=True)
    importer.add_argument("--destination", default="data/raw")
    importer.add_argument("--sha256")
    importer.add_argument("--member", help="Exact ZIP member name when an archive is ambiguous")
    importer.add_argument("--retrieval-date", help="YYYY-MM-DD; omit when unknown")
    importer.add_argument("--registry")
    importer.set_defaults(func=cmd_import)

    validate = data_sub.add_parser("validate", help="Parse and QC a JCAMP EI-MS file")
    validate.add_argument("--input", required=True)
    validate.add_argument("--source", required=True)
    validate.add_argument("--output")
    validate.add_argument("--registry")
    validate.set_defaults(func=cmd_validate)

    overlap = data_sub.add_parser(
        "audit-enfsi-overlap",
        help=(
            "Audit locally supplied ENFSI JDX spectra against SWGDRUG for spectral leakage; "
            "aggregate output contains no source spectra or compound names"
        ),
    )
    overlap.add_argument("--enfsi", required=True, help="ENFSI JDX/DX file or single-member ZIP")
    overlap.add_argument("--swgdrug", required=True, help="SWGDRUG HPJ/JDX spectrum file")
    overlap.add_argument(
        "--entry-list", help="Optional ENFSI entry-list XLSX for joint metadata audit"
    )
    overlap.add_argument(
        "--swgdrug-dataset",
        help="Normalized SWGDRUG JSONL; required with --entry-list",
    )
    overlap.add_argument("--near-identical-threshold", type=float, default=0.9999)
    overlap.add_argument("--output", help="Aggregate JSON report")
    overlap.add_argument(
        "--private-output",
        help=(
            "Optional local-only JSONL of record IDs/scores. Keep outside the open repository "
            "because ENFSI redistribution is restricted."
        ),
    )
    overlap.add_argument("--registry")
    overlap.set_defaults(func=cmd_data_audit_enfsi_overlap)

    external = data_sub.add_parser(
        "audit-enfsi-external",
        help=(
            "Construct conservative aggregate ENFSI external-test candidates after spectral, "
            "identity, and scaffold leakage filtering"
        ),
    )
    external.add_argument("--enfsi", required=True)
    external.add_argument("--swgdrug", required=True)
    external.add_argument("--entry-list", required=True)
    external.add_argument("--swgdrug-dataset", required=True)
    external.add_argument("--swgdrug-compounds", required=True)
    external.add_argument("--chebi-sdf", required=True)
    external.add_argument(
        "--supplemental-structures",
        default=str(_default_external_structures_path()),
        help=(
            "Public structure-evidence CSV used only for exact name/alias + formula "
            "fallback when ChEBI does not resolve a record"
        ),
    )
    external.add_argument(
        "--taxonomy-evidence",
        default=str(_default_evidence_path()),
        help="Public exact-match opioid/non-opioid evidence used for aggregate candidate labels",
    )
    external.add_argument("--near-identical-threshold", type=float, default=0.9999)
    external.add_argument("--output", required=True, help="Aggregate-only JSON report")
    external.add_argument(
        "--private-output",
        help=(
            "Optional local-only row manifest. Keep outside the open repository because "
            "ENFSI redistribution is restricted."
        ),
    )
    external.add_argument("--registry")
    external.set_defaults(func=cmd_data_audit_enfsi_external)

    builder = data_sub.add_parser("build", help="Build a versioned normalized dataset")
    builder.add_argument("--input", required=True)
    builder.add_argument("--source", required=True)
    builder.add_argument("--output", required=True)
    builder.add_argument("--dataset-version", default="opioidbench-eims-v0.1.0")
    builder.add_argument("--taxonomy", default="metadata/opioid_taxonomy.csv")
    builder.add_argument("--source-archive")
    builder.add_argument(
        "--compound-list",
        help="Optional SWGDRUG companion XLSX used for full names and version-new flags",
    )
    builder.add_argument(
        "--exact-mass-list",
        help=(
            "Optional SWGDRUG formula/exact-mass XLSX; values are retained as "
            "source-provided, formula-derived metadata"
        ),
    )
    builder.add_argument(
        "--agilent-archive",
        help=(
            "Optional SWGDRUG Agilent-format ZIP from the same release; only the embedded "
            "HEADER metadata is used to restore full source comments"
        ),
    )
    builder.add_argument(
        "--nist-archive",
        help=(
            "Optional SWGDRUG NIST-format ZIP from the same release; library-linked "
            "molfiles are mapped only after conservative formula/identity checks"
        ),
    )
    builder.add_argument("--registry")
    builder.set_defaults(func=cmd_build)

    metadata = sub.add_parser("metadata", help="Chemical metadata curation and enrichment")
    metadata_sub = metadata.add_subparsers(dest="metadata_command", required=True)
    enrich = metadata_sub.add_parser(
        "enrich",
        help="Resolve cached/PubChem structures by valid CAS and calculate RDKit descriptors",
    )
    enrich.add_argument("--input", required=True, help="Compound inventory CSV")
    enrich.add_argument("--output", required=True)
    enrich.add_argument("--cache", default="data/external/pubchem_cache")
    enrich.add_argument("--allow-network", action="store_true")
    enrich.add_argument("--max-records", type=int)
    enrich.set_defaults(func=cmd_metadata_enrich)

    resolve = metadata_sub.add_parser(
        "resolve",
        help=(
            "Resolve chemical identities with a local ChEBI bulk SDF first and PubChem "
            "CAS fallback second; unresolved rows are preserved for review"
        ),
    )
    resolve.add_argument("--input", required=True, help="Compound inventory CSV")
    resolve.add_argument("--output", required=True)
    resolve.add_argument(
        "--chebi-sdf",
        help="Local ChEBI FULL/LITE .sdf or .sdf.gz; no automatic bulk download is required",
    )
    resolve.add_argument("--cache", default="data/external/pubchem_cache")
    resolve.add_argument("--allow-network", action="store_true")
    resolve.add_argument("--max-records", type=int)
    resolve.set_defaults(func=cmd_metadata_resolve)

    taxonomy = metadata_sub.add_parser(
        "curate-taxonomy",
        help=(
            "Create curated taxonomy rows by exact name/CAS matching only; "
            "unmatched/conflicting compounds go to a local review queue"
        ),
    )
    taxonomy.add_argument("--input", required=True, help="Compound inventory CSV")
    taxonomy.add_argument(
        "--evidence", help="Evidence CSV; defaults to repository/package authoritative seed"
    )
    taxonomy.add_argument(
        "--exclusions",
        help="Exact target-exclusion CSV; defaults to repository/package exclusion seed",
    )
    taxonomy.add_argument("--output", default="metadata/opioid_taxonomy.csv")
    taxonomy.add_argument(
        "--review-output",
        default="results/private/opioid_taxonomy_review.csv",
    )
    taxonomy.set_defaults(func=cmd_metadata_curate_taxonomy)

    priority = metadata_sub.add_parser(
        "prioritize-taxonomy-review",
        help=(
            "Rank unresolved structured compounds by Morgan similarity to curated opioids; "
            "similarity never assigns a label"
        ),
    )
    priority.add_argument("--input", required=True, help="Resolved compound inventory CSV")
    priority.add_argument("--taxonomy", default="metadata/opioid_taxonomy.csv")
    priority.add_argument("--output", default="results/private/opioid_taxonomy_priority_review.csv")
    priority.add_argument("--threshold", type=float, default=0.75)
    priority.set_defaults(func=cmd_metadata_prioritize_review)

    enfsi = metadata_sub.add_parser(
        "audit-enfsi",
        help=(
            "Aggregate-only cross-check against a locally supplied ENFSI-DWG MS entry list; "
            "does not export restricted row-level source data"
        ),
    )
    enfsi.add_argument("--entry-list", required=True)
    enfsi.add_argument("--dataset", required=True, help="Normalized OpioidBench JSONL dataset")
    enfsi.add_argument("--output")
    enfsi.set_defaults(func=cmd_metadata_audit_enfsi)

    benchmark = sub.add_parser(
        "benchmark",
        help="Leakage-aware split generation and pilot baselines",
    )
    benchmark_sub = benchmark.add_subparsers(dest="benchmark_command", required=True)

    split = benchmark_sub.add_parser("split", help="Generate frozen cross-validation assignments")
    split.add_argument("--dataset", required=True, help="Normalized records JSONL")
    split.add_argument("--compounds", required=True, help="Resolved compound inventory CSV")
    split.add_argument("--taxonomy", default="metadata/opioid_taxonomy.csv")
    split.add_argument(
        "--protocol",
        choices=["random", "compound", "scaffold", "analogue"],
        required=True,
    )
    split.add_argument("--folds", type=int, default=5)
    split.add_argument("--seed", type=int, default=20260929)
    split.add_argument("--near-identical-threshold", type=float, default=0.9999)
    split.add_argument("--output", required=True)
    split.set_defaults(func=cmd_benchmark_split)

    baselines = benchmark_sub.add_parser("baselines", help="Run fixed pilot EI-MS baselines")
    baselines.add_argument("--dataset", required=True, help="Normalized records JSONL")
    baselines.add_argument("--compounds", required=True, help="Resolved compound inventory CSV")
    baselines.add_argument("--taxonomy", default="metadata/opioid_taxonomy.csv")
    baselines.add_argument("--assignments", required=True)
    baselines.add_argument(
        "--protocol",
        choices=["random", "compound", "scaffold", "analogue"],
        required=True,
    )
    baselines.add_argument("--mz-max", type=int, default=700)
    baselines.add_argument("--seed", type=int, default=20260929)
    baselines.add_argument("--output", required=True)
    baselines.set_defaults(func=cmd_benchmark_baselines)

    demo = sub.add_parser("demo", help="Run the synthetic first-milestone pipeline")
    demo.add_argument("--fixture", default="tests/fixtures/synthetic_eims.jdx")
    demo.add_argument("--output", default="results/demo")
    demo.add_argument("--registry")
    demo.set_defaults(func=cmd_demo)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    _configure_logging(args.verbose)
    try:
        return int(args.func(args))
    except (
        DataAccessPolicyError,
        SpectrumParseError,
        FileNotFoundError,
        KeyError,
        ValueError,
        RuntimeError,
    ) as exc:
        LOGGER.error("%s", exc)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
