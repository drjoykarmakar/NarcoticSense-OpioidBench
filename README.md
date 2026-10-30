# NarcoticSense-OpioidBench

**An open, structure-aware, uncertainty-aware benchmark for AI detection of opioids and emerging synthetic opioids from analytical spectra.**

> Status: Python package v0.4.0 / review-closed pilot v0.4.0. EI-MS ingestion/QC, chemical curation, leakage-aware pilot splits, and fixed classical pilot baselines are implemented. Reported metrics are preliminary, restricted to the evidence-curated subset, and not externally validated.

NarcoticSense-OpioidBench is the repository/benchmark, `opioidbench` is the Python package, and OpioidSpec-AI is the modeling framework planned for later milestones.

## Scientific motivation

Conventional spectral identification assumes that an unknown compound has an appropriate reference spectrum in the library. This project investigates the harder setting in which emerging opioid analogues may be absent from the training/reference library, while keeping exact identification, family recognition, open-set recognition, OOD detection, and uncertainty as separate tasks.

## Current milestone

Implemented now:

- reproducible Python package foundation and CI;
- source registry and source-adapter interface;
- conservative local import/download policy with SHA-256 provenance;
- ZIP-aware import for locally obtained JCAMP/HPJ source archives;
- EI-MS JCAMP/HPJ parser validated against the SWGDRUG 3.14 concatenated record layout;
- explicit raw -> validated -> normalized record layers;
- automated QC with machine-readable and human-readable reports;
- conservative CAS/name identity normalization without fuzzy merging;
- compound inventory and versioned dataset manifests;
- optional RDKit structure descriptors and explicit Morgan fingerprint configuration;
- conservative extraction/alignment of SWGDRUG NIST-format molecular structures with formula and isotope checks;
- streaming ChEBI bulk-SDF exact matching plus cached PubChem fallback for conservative chemical resolution;
- evidence-backed opioid/non-opioid-comparator taxonomy curation from authoritative WHO/UNODC public sources with no keyword auto-labeling;
- deterministic compound-, exact Murcko-scaffold-, and opioid-analogue-lineage-held-out 5-fold pilot splits with automated leakage auditing;
- fixed EI-MS pilot baselines: spectral 1-NN, logistic regression, linear SVM, and random forest;
- cached PubChem PUG REST enrichment, offline by default and rate-limited when enabled;
- SWGDRUG Agilent-header cross-format alignment to recover full source/contributor comments;
- aggregate-only ENFSI-DWG 2025.11 metadata, spectrum-overlap, identity, scaffold, and external-candidate auditing with restricted-source safeguards;
- synthetic, legally distributable fixtures and tests.

Not implemented yet: complete structure coverage, complete opioid/non-opioid taxonomy coverage, frozen final benchmark splits over the full intended label space, deep learning, OOD/open-set models, calibration, a sufficiently labeled ENFSI external test cohort, or manuscript-level final performance claims.

## Installation

```bash
python -m pip install -e '.[dev]'
# For RDKit-backed structure work:
python -m pip install -e '.[chem,dev]'
# Minimal extras for ENFSI spreadsheet + spectral overlap auditing:
python -m pip install -e '.[excel,spectral]'
# Leakage-aware split generation and classical pilot baselines:
python -m pip install -e '.[models]'
```

## Quick start

```bash
opioidbench data sources

# Import a legally obtained SWGDRUG archive. Raw data remain git-ignored.
opioidbench data import \
  --source swgdrug_eims \
  --input /path/to/JCAMP_01222025.zip \
  --destination data/raw/swgdrug_3.14 \
  --retrieval-date YYYY-MM-DD

# Build a versioned normalized dataset and QC audit.
opioidbench data build \
  --input 'data/raw/swgdrug_3.14/SWGDRUG 314.HPJ' \
  --source swgdrug_eims \
  --source-archive /path/to/JCAMP_01222025.zip \
  --output data/processed/opioidbench-eims-v0.1.0 \
  --dataset-version opioidbench-eims-v0.1.0

# Metadata enrichment is cache-only unless --allow-network is explicit.
opioidbench metadata enrich \
  --input data/processed/opioidbench-eims-v0.1.0/compounds.csv \
  --output data/processed/opioidbench-eims-v0.1.0/enrichment \
  --cache data/external/pubchem_cache

# Preferred chemical resolver: local ChEBI bulk SDF first, then PubChem fallback.
opioidbench metadata resolve \
  --input data/processed/opioidbench-eims-v0.1.0/compounds.csv \
  --chebi-sdf /path/to/chebi.sdf.gz \
  --output data/processed/opioidbench-eims-v0.1.0/resolution \
  --cache data/external/pubchem_cache

# Curate only exact authoritative evidence matches; unresolved rows go to a local review queue.
opioidbench metadata curate-taxonomy \
  --input data/processed/opioidbench-eims-v0.1.0/resolution/resolved_compounds.csv \
  --evidence metadata/opioid_evidence.csv \
  --exclusions metadata/taxonomy_exclusions.csv \
  --output metadata/opioid_taxonomy.csv \
  --review-output results/private/opioid_taxonomy_review.csv

# Rank unresolved structured compounds for human review; this never assigns a label.
opioidbench metadata prioritize-taxonomy-review \
  --input data/processed/opioidbench-eims-v0.1.0/resolution/resolved_compounds.csv \
  --taxonomy metadata/opioid_taxonomy.csv \
  --threshold 0.75 \
  --output results/private/opioid_taxonomy_priority_review.csv

# Optional restricted-source metadata cross-check; no ENFSI row data are exported.
opioidbench metadata audit-enfsi \
  --entry-list /path/to/2025ENFSI._DWG_MS_EntryList_2025_11.xlsx \
  --dataset data/processed/opioidbench-eims-v0.1.0/normalized_records.jsonl \
  --output results/private/enfsi_crosscheck.json

# Spectrum-level leakage audit against a locally supplied ENFSI JDX archive.
opioidbench data audit-enfsi-overlap \
  --enfsi /path/to/2025ENFSI-JDX.zip \
  --swgdrug 'data/raw/swgdrug_3.14/SWGDRUG 314.HPJ' \
  --entry-list /path/to/2025ENFSI._DWG_MS_EntryList_2025_11.xlsx \
  --swgdrug-dataset data/processed/opioidbench-eims-v0.1.0/normalized_records.jsonl \
  --output results/private/enfsi_spectral_overlap.json

# Freeze a compound-held-out pilot split. Use --protocol scaffold for the harder split.
opioidbench benchmark split \
  --dataset data/processed/opioidbench-eims-v0.1.0/normalized_records.jsonl \
  --compounds data/processed/opioidbench-eims-v0.1.0/resolution/resolved_compounds.csv \
  --taxonomy metadata/opioid_taxonomy.csv \
  --protocol compound \
  --output results/private/pilot_benchmark/compound_cv5

# Run the fixed pilot baselines on the frozen assignments.
opioidbench benchmark baselines \
  --dataset data/processed/opioidbench-eims-v0.1.0/normalized_records.jsonl \
  --compounds data/processed/opioidbench-eims-v0.1.0/resolution/resolved_compounds.csv \
  --taxonomy metadata/opioid_taxonomy.csv \
  --assignments results/private/pilot_benchmark/compound_cv5/compound_cv5_assignments.csv \
  --protocol compound \
  --output results/private/pilot_benchmark/compound_cv5/baselines

pytest
ruff check .
```

## SWGDRUG 3.14 local-build audit

Using the supplied SWGDRUG 3.14 JCAMP archive, the current parser/build pipeline produced:

- 3,826 spectra parsed and 3,826 passing peak/schema QC;
- 3,801 conservative compound identities;
- 1,502 records carrying the source placeholder CAS `000000-00-0`;
- 1 non-placeholder CAS value failing checksum validation;
- 2 exact duplicate-spectrum groups;
- 24 repeated conservative chemical-identity groups;
- 19 identity groups with conflicting source names;
- companion-list alignment, when supplied: 3,805 exact name matches, 19 verified
  fixed-width HPJ truncations, 2 preserved source disagreements, and 229 entries marked new in v3.14;
- formula/exact-mass workbook alignment, when supplied: all 3,826 formulas and names align to the
  source release, all 3,826 records receive a source-provided monoisotopic mass, 2,843 rows also
  carry the workbook's M+proton value, and 983 rows leave that optional value blank;
- 2 records (mercuric bromide and mercuric chloride) where rounding the source exact-mass value
  does not reproduce the HPJ nominal molecular weight; these remain warnings rather than corrections.
- Agilent-format header alignment, when supplied: all 3,826 records align; 1,779 source
  comments are identical to HPJ and 2,047 HPJ comments are verified prefixes of longer Agilent
  comments, allowing full contributor/provenance text to be retained without overwriting the raw HPJ.
- NIST-format structure alignment, when supplied: 3,297 library-linked structure blocks are found,
  3,296 parse with RDKit, and 3,207 spectra receive a conservatively accepted structure covering
  3,199 of 3,801 compound identities; 3,088 mappings have exact source/structure formula agreement
  and 119 are explicitly marked as unlabeled parent structures for deuterated isotopologues. The
  accepted structures yield 747 spectrum-level Bemis-Murcko scaffold IDs. At compound level,
  3,197 identities have one unambiguous canonical structure and 745 scaffold IDs; two source-CAS
  identity collisions contain conflicting structures and are excluded from compound/scaffold use.
  One spectrum with two distinct candidate structures is quarantined, 619 spectra remain without an
  accepted NIST structure, and four source-exact-mass/structure-mass differences are retained as
  warnings rather than used to override formula-consistent structures.
- ChEBI FULL bulk resolution, when supplied: 199,606 ChEBI entries are streamed locally; 938 of
  3,801 SWGDRUG compound identities resolve uniquely after exact InChIKey/CAS/name plus formula
  consistency checks. Of those, 844 corroborate an existing NIST structure and 94 fill previously
  missing structures. Seven rows remain ambiguous and 264 are rejected for consistency/structure
  conflicts. Final unambiguous structure coverage is 3,291/3,801 compounds (86.6%), leaving 510
  structureless identities. The ChEBI additions contribute 47 scaffold groups not present in the
  NIST-only compound inventory, raising the total exact scaffold-group count from 745 to 792 while
  preserving every pre-existing NIST scaffold identifier.

A locally supplied ENFSI-DWG 2025.11 entry list and JDX library were also checked without exporting
restricted row-level source data. The JDX contains 2,050/2,050 valid spectra and aligns with the
entry list by library-entry number and molecular formula for all 2,050 rows. It contains 18 exact
duplicate-spectrum pairs (36 records), 931 CAS `0` placeholders, and 613 additional CAS strings that
fail conservative CAS validation. Metadata alone maps 938
entry rows with high confidence to 855 SWGDRUG compounds; nine are ambiguous and 1,103 are unmatched.
A scale-invariant raw-intensity cosine audit finds zero byte/exact peak-list matches because ENFSI and
SWGDRUG use different intensity scales, but 103 ENFSI spectra have identical m/z support to a SWGDRUG
record and 361 have a best SWGDRUG cosine similarity >= 0.9999. Of those 361 near-identical leakage
risk flags, 88 occur among entries that the conservative metadata matcher left unmatched. Therefore
metadata-only deduplication is insufficient. The 1,689 records below the 0.9999 screen are only
*spectrum-disjoint candidates*; they are not yet proven compound-disjoint.

These are data-audit findings; source conflicts are not silently corrected. The current
authoritative evidence seed curates 147 opioid-positive compound identities and 121 source-defined
non-opioid comparators. Another 222 records are explicitly excluded from the biological target:
217 isotope-labelled/TMS analytical variants plus five exact target-policy exclusions. The previous
51-compound high-priority review queue is closed with 48 evidence-backed opioid decisions and three
N-benzyl target exclusions; recomputing the similarity frontier exposed eight additional close
neighbors, all of which were resolved with exact evidence. The remaining 3,311 identities are
unknown/pending, not negative, and the Morgan Tanimoto >=0.75 priority queue is now empty. This
supports controlled pilot experiments on the curated subset, but not a comprehensive
opioid-vs-non-opioid benchmark.


## Pilot leakage-aware baselines

A deterministic 5-fold pilot benchmark is frozen with seed `20260929` on one evidence snapshot.
The compound-held-out protocol uses 270 curated spectra (149 opioid, 121 comparator) and has zero
compound overlap, zero exact-spectrum overlap, and zero cross-fold raw-intensity cosine pairs >=
0.9999. The exact Murcko-scaffold-held-out protocol uses 257 structured curated spectra (139 opioid,
118 comparator) across 95 scaffold groups and has zero compound/scaffold overlap and zero >=0.9999
cross-fold spectral pairs.

The stricter analogue-lineage-held-out protocol uses those same 257 structured spectra but unions
exact scaffold groups with 17 broad curated opioid chemical lineages. Every fold has zero compound,
scaffold, opioid-lineage, exact-spectrum, and >=0.9999 spectral overlap. This grouping is an
evaluation constraint, not a new pharmacological label.

With nominal integer m/z bins through 700, square-root intensity transform, and L2 normalization,
the fixed random-forest baseline reaches mean balanced accuracy 0.938 and AUROC 0.987 on the
compound-held-out pilot, 0.850 and 0.878 on scaffold-held-out CV, and 0.616 and 0.651 on the stricter
analogue-lineage-held-out CV. The large drop under analogue holdout remains the central pilot result:
novel-chemistry performance cannot be inferred from the easier compound split.

A conservative ENFSI feasibility audit starts from 1,015 records that are both metadata-unmatched to
SWGDRUG and below the 0.9999 spectral-overlap screen. ChEBI plus three exact-name-and-formula public
structure fallbacks resolve 69 structures; exact InChIKey/scaffold de-leakage and internal duplicate
filtering leave 23 unique identity/scaffold candidates. Only two currently have source-backed opioid
labels, 20 have source-backed comparator labels, and one remains pending, so no external classifier
metric is reported yet.

These are preliminary subset/feasibility results, not final benchmark claims. See
`docs/benchmark_protocol.md` and `benchmarks/pilot_v0.4.0.md`.

## Data sources

The initial source is the SWGDRUG EI-MS library. The registry records provenance, observed checksums, and access policy. Because an explicit redistribution license has not been established, the repository does not include SWGDRUG raw spectra or processed spectrum payloads and does not automatically download them. Legally obtained local files can be imported with provenance and checksums.

See `docs/data_sources.md` and `metadata/source_registry.yaml`.

## Results

Aggregate pilot results and the external-test feasibility audit are reported in `benchmarks/pilot_v0.4.0.md`. Row-level real-source split
assignments and predictions remain local/git-ignored. Synthetic fixtures remain demonstrations of
plumbing only, and the pilot metrics are not final classification or external-detection claims.

## Reproducibility

`make reproduce` runs tests, lint checks, and the synthetic end-to-end demonstration without external research data. Real-source builds require a legally obtained local archive and remain outside Git.

## Disclaimer

Research-use computational screening framework. It is not independently validated for clinical, law-enforcement, evidentiary, diagnostic, or confirmatory analytical decision-making.

## License

Source code: BSD-3-Clause. External datasets retain their own terms and are not relicensed by this repository.


### SWGDRUG companion compound list

For SWGDRUG v3.14, pass the source-distributed companion workbook when available:

```bash
opioidbench data build \
  --input "data/raw/swgdrug_3.14/SWGDRUG 314.HPJ" \
  --source swgdrug_eims \
  --compound-list "/path/to/Compound List 01222025.xlsx" \
  --exact-mass-list "/path/to/SWGDrug 3.14 Excel list of names and formula.xlsx" \
  --agilent-archive "/path/to/AGILENT_01222025.zip" \
  --nist-archive "/path/to/NIST_01222025.zip" \
  --output data/processed/opioidbench-eims-v0.1.0
```

The build preserves both source names, uses the workbook only to resolve verified fixed-width
HPJ truncations, and records rows marked as new in source version 3.14.

### SWGDRUG formula/exact-mass companion

The source-distributed formula/exact-mass workbook may also be supplied with
`--exact-mass-list`. OpioidBench requires row/name/formula alignment before attaching its values.
The workbook itself states that the monoisotopic masses were calculated from the library formulas
with Agilent PCDL and are unverified; they are therefore stored as source-provided derived metadata,
not as independent measurements or curated structure evidence. The accompanying Agilent `.cdb`
file is recorded by checksum but is not parsed by the open-source pipeline.

### SWGDRUG Agilent-format provenance companion

`--agilent-archive` accepts the source-distributed SWGDRUG Agilent-format ZIP from the same
release. OpioidBench does **not** parse the proprietary vendor spectral binary. It reads only
printable metadata from the archive `HEADER`, validates release-order name/formula alignment
against the JCAMP build, and accepts a longer source comment only when the HPJ comment is an
exact prefix. This recovers contributor/provenance text while preserving the raw HPJ value.
The Agilent archive is another representation of SWGDRUG 3.14 and must not be counted as an
independent validation dataset.

### ENFSI-DWG external source

ENFSI-DWG version 2025.11 is registered as a separate candidate EI-MS source. Its public data
usage agreement restricts sharing/public release and third-party redistribution without prior
approval, so OpioidBench keeps ENFSI source files and row-level derivatives outside the open
repository. `metadata audit-enfsi` performs an aggregate metadata cross-check; `data
audit-enfsi-overlap` reads a local JDX/DX file or single-member ZIP and emits an aggregate spectral
overlap report. An optional local-only match manifest contains record IDs and similarity scores but
no names or peak lists and must remain outside the public repository. Because ENFSI/DWG is also a
documented contributor to SWGDRUG, even the post-screen subset must not be called an independent
compound-held-out test set until chemical identities/structures are resolved and leakage auditing is
complete.
