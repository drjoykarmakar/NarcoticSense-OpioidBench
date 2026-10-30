# Changelog

## 0.4.0 — review-closed pilot (first GitHub-ready source snapshot)

This version consolidates the pre-publication development milestones below. It is a pilot
research-code release, not a claim of externally validated classifier performance.

### Initial repository foundation

#### Added
- First milestone: repository foundation, source registry, safe local data import, EI-MS JCAMP parsing, schemas, QC, CLI, synthetic fixtures, tests, and CI.
- Real SWGDRUG 3.14 HPJ support for concatenated records without `##END=` markers.
- ZIP-aware local import with archive/member SHA-256 provenance.
- Versioned real-dataset build outputs, compound inventory, machine-readable audit, and human-readable QC summary.
- Conservative CAS/name chemical identity normalization and source-conflict detection.
- Optional RDKit structure descriptors and explicit Morgan fingerprint configuration.
- Human-curated opioid taxonomy schema and loader; no keyword auto-labeling.
- Cache-first PubChem PUG REST enrichment with explicit network opt-in, rate limiting, retry/backoff, and ambiguity rejection.

- Strict SWGDRUG companion-XLSX alignment for full source names and v3.14 new-entry flags.
- Separate preservation of HPJ and companion compound names, with repair limited to the
  observed 65-character HPJ truncation pattern and other source disagreements retained for review.

- Strict SWGDRUG formula/exact-mass workbook alignment with per-record source monoisotopic mass,
  optional source M+proton values, checksum provenance, and source-value consistency warnings.
- Registry provenance for the supplied `EXACTMASS_01222025.zip`, its XLSX member, and Agilent PCDL
  `.cdb` member; the proprietary CDB is recorded but not parsed.

### source cross-check extensions

- Added conservative extraction of full SWGDRUG source comments from the same-release Agilent
  `HEADER` metadata, with strict JCAMP name/formula/comment-prefix validation.
- Added an ENFSI-DWG 2025.11 registry entry and aggregate-only cross-source audit command that does
  not export restricted ENFSI row-level data.
- Recorded supplied Agilent archive/header and ENFSI entry-list SHA-256 checksums.

### ENFSI spectrum-overlap audit

- Added restricted-source ENFSI JDX/ZIP spectrum ingestion for local auditing without persisting the
  supplied JDX into repository artifacts.
- Added scale-invariant cross-library cosine overlap screening, exact peak-support checks, threshold
  sensitivity, and an optional local-only ID/score manifest.
- Validated the supplied ENFSI 2025.11 JDX: 2,050/2,050 spectra valid; 103 records share identical
  m/z support with SWGDRUG 3.14 and 361 have best cosine >= 0.9999.
- Joint metadata/spectrum auditing shows that 88 of the >=0.9999 spectral overlap flags are missed by
  the conservative metadata matcher, so metadata-only external-set deduplication is insufficient.
- Recorded the supplied ENFSI JDX archive/member SHA-256 checksums while retaining restricted source
  files outside Git.

### SWGDRUG NIST structure extraction

- Added conservative parsing of library-linked molfiles from the same-release NIST `USRSTRUC.DB`,
  including the older pre-version-tag molfile variant present in the archive.
- Added formula-aware mapping by current library ID with exact-name/CAS fallback, explicit
  `isotopologue_parent` handling for deuterated records, and quarantine of competing structures.
- On the supplied SWGDRUG 3.14 files, 3,207/3,826 spectra receive accepted structures covering
  3,199/3,801 compound identities and 747 exact Bemis-Murcko scaffold IDs.
- Added local NIST structure review output and archive/member SHA-256 provenance; raw/vendor
  structure payloads remain outside Git.
- Recorded the supplied Shimadzu archive and SWGDRUG changelog checksums as provenance artifacts.

### ChEBI resolution and authoritative taxonomy seed

- Added streaming local ChEBI FULL/LITE SDF resolution with exact InChIKey/CAS/name keys,
  formula/exact-mass consistency checks, and ambiguity/conflict quarantine.
- Added ChEBI-first/PubChem-second compound resolution outputs and a local chemical review queue;
  network access remains explicit and optional.
- Added a public opioid evidence seed sourced from WHO ATC/DDD opioid groups and WHO ECDD reports,
  with pharmacological evidence separated from alternate-name/CAS identity corroboration.
- At this initial taxonomy-seed milestone, curated 37 exact opioid-positive SWGDRUG compound
  identities across nine source-defined groups; 3,764 remained pending and no non-opioid comparator
  seed had yet been added.
- Kept pending SWGDRUG-derived rows out of the tracked taxonomy file to avoid republishing the source
  compound inventory while redistribution rights remain unresolved.
- Clarified structure coverage: 3,199 compound IDs have structure evidence, 3,197 have one
  unambiguous compound-level structure, and two source-CAS collisions contain conflicting structures;
  745 scaffold IDs remain at the unambiguous compound level.


### real ChEBI bulk resolution

- Ran the conservative resolver against the supplied ChEBI FULL SDF snapshot: 199,606 entries
  scanned and 938/3,801 SWGDRUG compound identities accepted as exact, consistency-checked matches.
- Added 94 structures beyond the NIST-only inventory, increasing unambiguous compound structure
  coverage from 3,197 to 3,291 (86.6%) and leaving 510 identities structureless.
- Added 47 new scaffold groups beyond the NIST-only inventory, increasing the exact compound-level
  scaffold-group count from 745 to 792.
- Quarantined seven ambiguous ChEBI matches and 264 consistency-conflict matches rather than
  overriding source structures or formulas.
- Preserved previously accepted source-linked NIST scaffold identifiers during external resolution;
  re-parsing existing SMILES is no longer allowed to silently rewrite scaffold hashes.
- Recorded the supplied ChEBI snapshot SHA-256 and CC BY 4.0 provenance in the source registry.

### source-backed non-opioid comparators and priority review

- Expanded the public evidence seed with WHO ATC/DDD 2026 comparator classes and WHO ECDD
  source-defined NPS categories while retaining exact CAS/name/alias matching only.
- Added additional WHO ECDD/fact-sheet opioid evidence plus six supplied ChEBI Release 255 records
  whose definitions explicitly state that the compound is a synthetic opioid.
- Current real-data curation yields 55 opioid positives, 118 source-defined non-opioid comparators,
  3,628 unknown/pending identities, and zero evidence conflicts.
- Added `metadata prioritize-taxonomy-review`, which ranks unresolved structured compounds by
  chirality-aware Morgan Tanimoto similarity to curated opioids without assigning labels. At the
  default 0.75 threshold, 89 unresolved compounds are prioritized for human review.
- Added regression tests for explicit non-opioid evidence, conflicting-evidence quarantine, and
  similarity-based review prioritization.

### leakage-aware pilot benchmark

- Added ten exact public WHO/UNODC evidence records from the high-priority review queue, increasing
  curated opioid positives from 55 to 66 while retaining 118 source-defined non-opioid comparators;
  3,617 identities remain unknown/pending and zero evidence conflicts are present.
- Added deterministic random/compound/scaffold cross-validation split generation with fold-level
  compound, exact Murcko scaffold, exact-spectrum, and >=0.9999 cross-fold cosine leakage audits.
- Froze 5-fold pilot compound-held-out and scaffold-held-out protocols at seed `20260929`. The
  compound protocol contains 184 curated spectra; the scaffold protocol contains 173 structured
  curated spectra across 86 exact scaffold groups. Both have zero exact-spectrum and >=0.9999
  cross-fold spectral overlap; the scaffold protocol has zero exact scaffold overlap in every fold.
- Added fixed spectral 1-NN, logistic regression, linear SVM, and random-forest EI-MS baselines with
  nominal m/z binning, square-root intensity transform, and L2 normalization.
- Added aggregate-only pilot reporting; real assignments/predictions remain ignored local outputs.

### analogue holdout and ENFSI chemical de-leakage

- Expanded exact evidence curation to 72 opioid positives and 120 source-defined non-opioid
  comparators; 3,609 SWGDRUG identities remain unknown/pending with zero evidence conflicts.
- Added deterministic opioid analogue-lineage-held-out CV that unions exact Murcko-scaffold groups
  with broad curated opioid chemical lineages; every fold has zero compound, scaffold, lineage,
  exact-spectrum, and >=0.9999 spectral overlap.
- Refreshed all three leakage-aware pilot protocols on the same evidence snapshot: 192 spectra for
  compound holdout and 181 structured spectra for scaffold/analogue holdout.
- Added exact-name-plus-formula public structure fallbacks for ENFSI records unresolved by ChEBI;
  supplemental structures are never selected by CAS alone.
- Added aggregate taxonomy labeling to the ENFSI external-candidate audit while keeping all
  row-level ENFSI derivatives private/restricted.
- The external feasibility cascade reduces 1,015 metadata-unmatched/spectrum-disjoint ENFSI records
  to 23 unique identity/scaffold candidates: two source-backed opioids, 12 source-backed
  non-opioid comparators, and nine pending. No external performance metric is reported from this
  undersized positive set.
- Added regression tests for supplemental structure guards, resource synchronization, and analogue
  component grouping.

### target-clean taxonomy and pilot v0.3

- Expanded exact public evidence curation to 88 opioid positives and 121 source-backed non-opioid
  comparators with zero authoritative label conflicts.
- Added `metadata/taxonomy_exclusions.csv` and package-synchronized target-exclusion handling. The
  current real build excludes 217 isotope-labelled/TMS analytical variants plus dextrorphan and
  6beta-naltrexol from the biological target instead of forcing them into either class.
- Reduced the structure-similarity priority review queue to 51 unresolved compounds; 3,373 source
  identities remain unknown/pending and are not treated as negative.
- Refreshed all three frozen five-fold protocols at seed `20260929`: 209 curated spectra for
  compound holdout and 198 structured spectra for scaffold/analogue-lineage holdout. All required
  leakage constraints remain zero for their respective protocols.
- Pilot v0.3 random-forest mean balanced accuracy/AUROC is 0.920/0.979 (compound), 0.742/0.848
  (scaffold), and 0.552/0.713 (analogue lineage). Spectral 1-NN has the highest analogue-held-out
  balanced accuracy among the four fixed baselines at 0.609.
- Refreshed the ENFSI aggregate feasibility audit against the rebuilt SWGDRUG dataset. The 23
  chemically de-leaked candidates now comprise two source-backed opioids, 20 source-backed
  non-opioid comparators, and one pending identity; no external performance metric is reported.

### review-closed taxonomy and pilot v0.4

- Closed the prior 51-compound structure-prioritized review queue: 48 exact evidence-backed opioid
  labels and three N-benzyl fentanyl-related target exclusions, with no similarity-derived labels.
- Recomputed the review frontier after each curation pass; eight newly exposed near neighbors were
  also resolved with exact evidence. The Morgan Tanimoto >=0.75 priority queue is now empty.
- Expanded the exact public evidence seed to 147 opioid positives and 121 source-backed non-opioid
  comparators. The target-clean layer now excludes 217 analytical variants plus five exact
  target-policy cases; 3,311 source identities remain pending and zero evidence conflicts are present.
- Refreshed five-fold v0.4 protocols at seed `20260929`: 270 curated spectra for compound holdout and
  257 structured spectra for scaffold/analogue-lineage holdout. Required leakage constraints remain
  zero for their respective protocols.
- Pilot v0.4 random-forest mean balanced accuracy/AUROC is 0.938/0.987 (compound), 0.850/0.878
  (scaffold), and 0.616/0.651 (analogue lineage).
- Refreshed ENFSI external-test feasibility against the final evidence snapshot. The chemically
  de-leaked pool remains 23 unique identity/scaffold candidates: two source-backed opioids, 20
  source-backed comparators, and one pending identity; no external performance metric is reported.
- Added regression tests for exact N-benzyl target exclusions and frontier evidence aliases.
