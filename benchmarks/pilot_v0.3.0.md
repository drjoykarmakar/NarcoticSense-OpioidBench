# Pilot Benchmark v0.3.0

## Scope

This pilot uses only exact evidence-curated SWGDRUG labels: 88 opioid compounds and 121 source-defined non-opioid comparators. A further 219 records are explicitly excluded from the biological target (217 isotope-labelled/TMS analytical variants plus two exact pharmacology-based target exclusions). The other 3,373 conservative compound identities remain unknown/pending and are excluded from supervised evaluation.

All splits use five folds and seed `20260929`. Row-level assignments and predictions remain local/git-ignored.

## Frozen protocols

| Protocol | Spectra | Opioid | Comparator | Exact scaffolds | Opioid lineages | Enforced cross-fold constraints |
|---|---:|---:|---:|---:|---:|---|
| Compound-held-out | 209 | 88 | 121 | 92 observed | 16 observed | compound, exact spectrum, cosine >=0.9999 |
| Scaffold-held-out | 198 | 80 | 118 | 92 | 15 | compound, scaffold, exact spectrum, cosine >=0.9999 |
| Analogue-lineage-held-out | 198 | 80 | 118 | 92 | 15 | compound, scaffold, opioid analogue lineage, exact spectrum, cosine >=0.9999 |

Eleven curated spectra without an accepted exact scaffold are excluded from the two structure-aware protocols. The analogue protocol is an operational evaluation grouping, not a new pharmacological label.

## Fixed baselines

Spectra are represented on nominal integer m/z bins through 700, square-root transformed, and L2-normalized. The models are spectral 1-NN, logistic regression, linear SVM, and a 500-tree random forest with the fixed settings documented in `docs/benchmark_protocol.md`.

| Protocol | Model | Balanced accuracy mean | Balanced accuracy SD | AUROC mean | AUROC SD | Average precision mean |
|---|---|---:|---:|---:|---:|---:|
| Compound-held-out | linear svm | 0.870 | 0.062 | 0.942 | 0.049 | 0.916 |
| Compound-held-out | logistic regression | 0.837 | 0.097 | 0.901 | 0.088 | 0.861 |
| Compound-held-out | random forest | 0.920 | 0.060 | 0.979 | 0.026 | 0.972 |
| Compound-held-out | spectral 1nn | 0.878 | 0.045 | 0.921 | 0.054 | 0.877 |
| Scaffold-held-out | linear svm | 0.675 | 0.189 | 0.711 | 0.252 | 0.690 |
| Scaffold-held-out | logistic regression | 0.651 | 0.160 | 0.653 | 0.301 | 0.618 |
| Scaffold-held-out | random forest | 0.742 | 0.172 | 0.848 | 0.181 | 0.814 |
| Scaffold-held-out | spectral 1nn | 0.701 | 0.147 | 0.684 | 0.316 | 0.681 |
| Analogue-lineage-held-out | linear svm | 0.548 | 0.127 | 0.566 | 0.258 | 0.544 |
| Analogue-lineage-held-out | logistic regression | 0.543 | 0.141 | 0.536 | 0.258 | 0.483 |
| Analogue-lineage-held-out | random forest | 0.552 | 0.062 | 0.713 | 0.207 | 0.678 |
| Analogue-lineage-held-out | spectral 1nn | 0.609 | 0.190 | 0.639 | 0.301 | 0.628 |

## Interpretation

The compound-held-out random forest is strong on this small curated subset (balanced accuracy 0.920; AUROC 0.979), but random-forest performance drops under exact scaffold holdout (0.742; 0.848) and further in balanced accuracy under the stricter analogue-lineage holdout (0.552; 0.713). Under analogue holdout, spectral 1-NN has the highest balanced accuracy among these four fixed baselines (0.609). The large novelty gap means performance on unseen opioid chemistry must not be inferred from the easier compound-held-out result.

## ENFSI external-test feasibility

ENFSI source spectra and row-level derivatives remain local/restricted. Starting from 1,015 ENFSI records that are both metadata-unmatched to SWGDRUG and below the 0.9999 spectral-overlap threshold, conservative ChEBI matching plus three exact-name-and-formula public structure fallbacks resolves 69 structures. Exact InChIKey overlap, exact Murcko-scaffold overlap, and ENFSI internal exact-duplicate filtering leave 23 unique identity/scaffold candidates.

Of those 23 chemically de-leaked candidates, 2 have source-backed opioid labels, 20 have source-backed non-opioid comparator labels, and 1 remains pending. This is not enough positive support for a defensible external classifier metric, so v0.3.0 reports feasibility counts only and deliberately does not report an external AUROC or balanced accuracy.

## Reproducibility and limits

Row-level real-source assignments, predictions, taxonomy review queues, and ENFSI candidate manifests stay local and git-ignored. Public artifacts contain aggregate metrics and hashes only. These results are pilot subset results and do not establish performance on all SWGDRUG compounds, on a broad unseen-opioid distribution, or for operational/confirmatory forensic use.
