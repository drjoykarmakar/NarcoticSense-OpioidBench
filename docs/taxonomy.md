# Opioid Taxonomy

`metadata/opioid_taxonomy.csv` is a human-reviewable curation table with these fields:

```text
compound_id
preferred_name
opioid_status
opioid_family
subfamily
evidence_source
curation_status
notes
```

Allowed curation statuses are `pending`, `curated`, `uncertain`, and `excluded`.

## Operational labels

The benchmark uses three label states plus an explicit target-exclusion state:

- `opioid`: an authoritative source explicitly identifies the substance as an opioid, opioid
  analgesic/anesthetic, fentanyl analogue in an opioid context, or source-defined novel synthetic
  opioid.
- `non-opioid`: an authoritative source explicitly places the substance in a distinct comparator
  class or otherwise supports its use as a non-opioid benchmark comparator. This is an operational
  comparator label, not a claim of zero affinity at every opioid receptor.
- unknown/pending: no exact authoritative evidence match has been established. Absence from an
  opioid list never creates a negative label.
- `excluded`: the record is inappropriate for the biological target rather than positive/negative.
  This includes isotope-labelled/TMS analytical variants and exact target-policy exclusions.

WHO ATC/DDD, WHO ECDD, UNODC/INCB public material, FDA/NCATS records, institutional forensic opioid
reference collections, and peer-reviewed pharmacology are used as traceable evidence. ChEBI/PubChem
are used for exact identity/structure corroboration where needed.

## Matching policy

Evidence is matched by exact normalized name/alias or validated CAS only. Fuzzy name matching,
substring rules, scaffold membership, and spectral or structural similarity do not assign labels.
Conflicting exact evidence is sent to the review queue as uncertain.

The tracked evidence seed is `metadata/opioid_evidence.csv`. Exact biological-target exclusions are
tracked separately in `metadata/taxonomy_exclusions.csv`. Their packaged copies must be
byte-identical to the repository copies. Unresolved review data remain local/git-ignored because
they contain source-derived SWGDRUG row-level metadata.

Analytical variants are excluded conservatively from the supervised biological target when a
resolved molecular formula/structure shows isotope labelling or when the source identity is an
explicit TMS derivatization product. Exclusion is not label propagation from the parent compound.

## Current local audit

Against the 3,801 conservative SWGDRUG compound identities, the current seed yields:

- 147 curated opioid positives;
- 121 curated source-defined non-opioid comparators;
- 217 analytical-variant exclusions;
- 5 exact target-policy exclusions;
- 3,311 unknown/pending identities;
- zero authoritative label conflicts.

The previous 51-compound high-priority queue was reviewed individually. Forty-eight received exact
evidence-backed opioid labels; three N-benzyl fentanyl-related edge cases were target-excluded
rather than inferred positive or negative. Recomputing the similarity frontier after those decisions
exposed eight additional near-neighbor compounds, and all eight were resolved with exact evidence.

Of the remaining unresolved identities, 2,873 have an accepted structure. Using the same Morgan
fingerprint configuration as the chemistry pipeline (radius 2, 2,048 bits, chirality enabled), a
Tanimoto threshold of 0.75 now prioritizes zero unresolved compounds against 137 structured curated
opioid references. An empty priority queue does not mean the remaining 3,311 compounds have been
classified; it means none meet this specific structure-similarity review threshold.

Run:

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

Similarity never creates a label. Remaining identities that do not meet an explicit exclusion or
evidence rule stay pending.
