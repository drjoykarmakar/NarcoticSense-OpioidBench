# Pilot Benchmark v0.4.0

## Scope

This pilot uses only exact evidence-curated SWGDRUG labels: 147 opioid compounds and 121
source-defined non-opioid comparators. A further 222 records are explicitly excluded from the
biological target (217 isotope-labelled/TMS analytical variants plus five exact target-policy
exclusions). The other 3,311 conservative compound identities remain unknown/pending and are
excluded from supervised evaluation.

The v0.4 curation closes the previous 51-compound priority queue: 48 of those compounds received
exact evidence-backed opioid labels and three N-benzyl fentanyl-related edge cases were
intentionally target-excluded. Recomputing structure similarity after those decisions exposed eight
additional near-neighbor compounds; all eight were resolved with exact evidence. The review queue at
Morgan Tanimoto >= 0.75 is therefore empty. This does not imply that all 3,311 pending compounds are
biologically classified.

All splits use five folds and seed `20260929`. Row-level assignments and predictions remain
local/git-ignored.

## Frozen protocols

| Protocol | Spectra | Opioid | Comparator | Exact scaffolds | Opioid lineages | Enforced cross-fold constraints |
|---|---:|---:|---:|---:|---:|---|
| Compound-held-out | 270 | 149 | 121 | 95 observed | 18 observed | compound, exact spectrum, cosine >=0.9999 |
| Scaffold-held-out | 257 | 139 | 118 | 95 | 17 | compound, scaffold, exact spectrum, cosine >=0.9999 |
| Analogue-lineage-held-out | 257 | 139 | 118 | 95 | 17 | compound, scaffold, opioid analogue lineage, exact spectrum, cosine >=0.9999 |

Thirteen curated spectra without an accepted exact scaffold are excluded from the two
structure-aware protocols. The analogue protocol is an operational evaluation grouping, not a new
pharmacological label.

## Fixed baselines

Spectra are represented on nominal integer m/z bins through 700, square-root transformed, and
L2-normalized. The models are spectral 1-NN, logistic regression, linear SVM, and a 500-tree random
forest with the fixed settings documented in `docs/benchmark_protocol.md`.

| Protocol | Model | Balanced accuracy mean | Balanced accuracy SD | AUROC mean | AUROC SD | Average precision mean |
|---|---|---:|---:|---:|---:|---:|
| Compound-held-out | linear svm | 0.900 | 0.037 | 0.971 | 0.016 | 0.977 |
| Compound-held-out | logistic regression | 0.865 | 0.045 | 0.942 | 0.026 | 0.953 |
| Compound-held-out | random forest | 0.938 | 0.035 | 0.987 | 0.013 | 0.990 |
| Compound-held-out | spectral 1nn | 0.888 | 0.040 | 0.970 | 0.024 | 0.977 |
| Scaffold-held-out | linear svm | 0.780 | 0.093 | 0.792 | 0.240 | 0.791 |
| Scaffold-held-out | logistic regression | 0.731 | 0.123 | 0.727 | 0.258 | 0.741 |
| Scaffold-held-out | random forest | 0.850 | 0.139 | 0.878 | 0.152 | 0.878 |
| Scaffold-held-out | spectral 1nn | 0.776 | 0.169 | 0.798 | 0.287 | 0.824 |
| Analogue-lineage-held-out | linear svm | 0.504 | 0.120 | 0.540 | 0.300 | 0.540 |
| Analogue-lineage-held-out | logistic regression | 0.524 | 0.141 | 0.540 | 0.300 | 0.547 |
| Analogue-lineage-held-out | random forest | 0.616 | 0.155 | 0.651 | 0.248 | 0.624 |
| Analogue-lineage-held-out | spectral 1nn | 0.519 | 0.138 | 0.571 | 0.279 | 0.547 |

## Interpretation

The compound-held-out random forest is strong on this curated subset (balanced accuracy 0.938;
AUROC 0.987). Exact scaffold holdout is harder (0.850; 0.878), and the strict
analogue-lineage-held-out protocol drops further (0.616; 0.651). The novelty gap remains the central
pilot result: performance on unseen opioid chemistry must not be inferred from compound-held-out
cross-validation.

## ENFSI external-test feasibility

ENFSI source spectra and row-level derivatives remain local/restricted. Starting from 1,015 ENFSI
records that are both metadata-unmatched to SWGDRUG and below the 0.9999 spectral-overlap threshold,
conservative ChEBI matching plus three exact-name-and-formula public structure fallbacks resolves 69
structures. Exact InChIKey overlap, exact Murcko-scaffold overlap, and ENFSI internal exact-duplicate
filtering leave 23 unique identity/scaffold candidates.

Of those 23 chemically de-leaked candidates, 2 have source-backed opioid labels, 20 have
source-backed non-opioid comparator labels, and 1 remains pending. Finishing the SWGDRUG priority
queue therefore does not solve the external-positive scarcity. This is not enough positive support
for a defensible external classifier metric, so v0.4.0 reports feasibility counts only and
deliberately does not report an external AUROC or balanced accuracy.

## Reproducibility and limits

Row-level real-source assignments, predictions, unresolved review rows, and ENFSI candidate
manifests stay local and git-ignored. Public artifacts contain aggregate metrics and hashes only.
These results are pilot subset results and do not establish performance on all SWGDRUG compounds, on
a broad unseen-opioid distribution, or for operational/confirmatory forensic use.
