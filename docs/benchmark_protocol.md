# Benchmark Protocol

## Status

OpioidBench includes a pilot leakage-aware evaluation protocol on the evidence-curated SWGDRUG
subset. It is intentionally not the final benchmark because 3,311 of 3,801 conservative compound
identities remain unknown/pending and the chemically de-leaked ENFSI candidate pool has only two
source-backed opioid positives.

## Pilot label set and target exclusions

Only taxonomy rows with `curation_status=curated` and `opioid_status` equal to `opioid` or
`non-opioid` enter supervised evaluation. Absence from an opioid source never creates a negative
label. Current counts are 147 opioid compounds and 121 source-defined comparator compounds.

Rows with `curation_status=excluded` are neither positive nor negative. The current exclusion layer
removes 217 isotope-labelled/TMS analytical variants plus five exact pharmacology/target-policy
exclusions. Three of those five are N-benzyl fentanyl-related edge cases for which exact
derivative-specific agonism is not established; they are excluded rather than inferred positive or
negative. This prevents analytical standards, derivatives, antagonists, and unresolved target-edge
cases from contaminating either class.

## Frozen cross-validation protocols

All pilot splits use five folds and seed `20260929`.

### Compound-held-out CV

All spectra from a `compound_id` are assigned to the same fold. Current local build: 270 spectra,
149 opioid and 121 comparator, representing 268 curated compounds. Every fold has zero train/test
compound overlap. Exact duplicate peak lists and raw-intensity cosine >= 0.9999 cross-fold pairs are
also audited and are zero in every fold.

### Exact Murcko-scaffold-held-out CV

Only curated compounds with an accepted, non-conflicting exact Murcko scaffold enter this protocol.
All compounds sharing a scaffold are assigned to one fold. Current local build: 257 spectra, 139
opioid and 118 comparator, spanning 95 exact scaffold groups. Thirteen curated spectra without an
exact scaffold are excluded. Every fold has zero train/test compound and scaffold overlap, zero
exact-spectrum overlap, and zero >=0.9999 cross-fold spectral pairs.

### Opioid analogue-lineage-held-out CV

This protocol starts from the same 257 structured curated spectra. It constructs connected split
groups by taking the union of exact Murcko-scaffold groups and broad curated opioid chemical
lineages. The lineages are operational evaluation groups, not automatically inferred opioid labels.
Their purpose is to prevent closely related curated opioid chemistry from appearing in both train
and test merely because the exact Murcko scaffolds differ.

The current build contains 17 observed opioid analogue lineages. Every fold has zero train/test
compound, scaffold, opioid-lineage, exact-spectrum, and >=0.9999 spectral overlap.

Acyclic/structureless records are not given artificial per-compound scaffold identifiers for these
structure-aware splits.

## Fixed pilot baselines

Spectra are binned to nominal integer m/z 0--700, intensities are square-root transformed, and each
vector is L2-normalized. Four fixed baselines are evaluated without fold-specific hyperparameter
tuning:

- spectral 1-NN: maximum cosine to each training class, decision = best opioid cosine minus best
  comparator cosine;
- logistic regression: C=1, balanced class weights;
- linear SVM: LinearSVC C=1, balanced class weights;
- random forest: 500 trees, square-root feature sampling, balanced-subsample class weights.

Metrics include accuracy, balanced accuracy, sensitivity, specificity, precision, F1, MCC, AUROC,
and average precision. Full aggregate pilot results are documented in `benchmarks/pilot_v0.4.0.md`.

## ENFSI external-test feasibility gate

Restricted ENFSI row-level data remain local. External candidates must be metadata-unmatched to
SWGDRUG, below the 0.9999 SWGDRUG spectral-overlap screen, structurally resolved, absent from
SWGDRUG by exact InChIKey and exact Murcko scaffold, and outside an ENFSI exact duplicate-spectrum
group. Public supplemental structures may fill a ChEBI-unresolved record only after exact
name/alias plus molecular-formula agreement; CAS alone cannot select a supplemental structure.

The current aggregate audit reduces 1,015 metadata-unmatched/spectrum-disjoint records to 23 unique
identity/scaffold candidates. Two have source-backed opioid labels, 20 have source-backed comparator
labels, and one remains pending. Therefore the project does not report an external AUROC, balanced
accuracy, or other classifier metric from ENFSI at this stage.

## Interpretation limits

The compound-held-out random-forest result is strong on the curated subset, but novelty constraints
remain material. Current random-forest mean balanced accuracy/AUROC are 0.938/0.987 for compound
holdout, 0.850/0.878 for scaffold holdout, and 0.616/0.651 for analogue-lineage holdout. The decline
under analogue holdout shows that performance on unseen opioid chemistry cannot be inferred from the
easier compound split.

These results characterize the evidence-curated subset only. They do not establish performance on
all SWGDRUG compounds, a broad unseen-opioid distribution, open-set/OOD recognition, or
operational/confirmatory forensic use.
