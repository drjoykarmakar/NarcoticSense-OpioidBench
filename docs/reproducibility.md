# Reproducibility

## Synthetic CI path

```bash
make setup
make reproduce
```

This runs linting, tests, and a synthetic EI-MS parse -> validate -> QC -> normalized-record demonstration. No external research data are required.

## SWGDRUG 3.14 local build

With a legally obtained `JCAMP_01222025.zip`:

```bash
opioidbench data import \
  --source swgdrug_eims \
  --input /path/to/JCAMP_01222025.zip \
  --destination data/raw/swgdrug_3.14 \
  --retrieval-date YYYY-MM-DD

opioidbench data build \
  --input 'data/raw/swgdrug_3.14/SWGDRUG 314.HPJ' \
  --source swgdrug_eims \
  --source-archive /path/to/JCAMP_01222025.zip \
  --output data/processed/opioidbench-eims-v0.1.0 \
  --dataset-version opioidbench-eims-v0.1.0
```

Expected observed checksums for the development archive are documented in `metadata/source_registry.yaml`. Generated raw/processed data are git-ignored.

The build writes `normalized_records.jsonl`, `data_qc.json`, `data_qc.md`, `compounds.csv`, `dataset_audit.json`, and `dataset_manifest.json`.

## Metadata enrichment

`opioidbench metadata enrich` is offline by default. Use `--allow-network` only when public API access is appropriate. Cache contents belong under `data/external/` and are git-ignored.

## Evidence curation and target exclusions

After chemical resolution, rebuild the public curation layer with exact evidence and target
exclusions:

```bash
opioidbench metadata curate-taxonomy \
  --input data/processed/opioidbench-eims-v0.1.0/resolution/resolved_compounds.csv \
  --evidence metadata/opioid_evidence.csv \
  --exclusions metadata/taxonomy_exclusions.csv \
  --output metadata/opioid_taxonomy.csv \
  --review-output results/private/opioid_taxonomy_review.csv

opioidbench metadata prioritize-taxonomy-review \
  --input data/processed/opioidbench-eims-v0.1.0/resolution/resolved_compounds.csv \
  --taxonomy metadata/opioid_taxonomy.csv \
  --threshold 0.75 \
  --output results/private/opioid_taxonomy_priority_review.csv
```

The current local snapshot yields 147 curated opioid compounds, 121 source-backed comparators, 222
target exclusions (217 analytical variants plus five exact target-policy exclusions), 3,311 pending
identities, and zero high-priority unresolved compounds at the 0.75 similarity threshold.

## Frozen pilot v0.4

Use seed `20260929`, five folds, and near-identical spectral threshold `0.9999`. Generate the
`compound`, `scaffold`, and `analogue` protocols with `opioidbench benchmark split`, then run the
fixed models with `opioidbench benchmark baselines`. Aggregate reference results are recorded in
`benchmarks/pilot_v0.4.0.md`; row-level assignments and predictions stay under `results/private/`.
