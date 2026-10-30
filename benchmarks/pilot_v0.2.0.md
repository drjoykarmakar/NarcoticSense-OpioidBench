# Pilot Benchmark v0.2.0

## Frozen label set

This pilot uses only exact evidence-curated SWGDRUG labels: 72 opioid compounds and 120 source-defined non-opioid comparators. The other 3,609 conservative compound identities remain unknown/pending and are excluded from supervised evaluation.

All cross-validation protocols use five folds and seed `20260929`. No hyperparameter tuning is performed within this pilot.

## Leakage-aware protocols

| Protocol | Spectra | Opioid | Comparator | Exact scaffolds | Opioid analogue lineages | Required zero-overlap constraints |
|---|---:|---:|---:|---:|---:|---|
| Compound-held-out | 192 | 72 | 120 | 91 observed | 15 observed | compound, exact spectrum, cosine >=0.9999 |
| Scaffold-held-out | 181 | 64 | 117 | 91 | 14 | compound, scaffold, exact spectrum, cosine >=0.9999 |
| Analogue-lineage-held-out | 181 | 64 | 117 | 91 | 14 | compound, scaffold, opioid analogue lineage, exact spectrum, cosine >=0.9999 |

The analogue protocol is an operational evaluation grouping, not a new pharmacological label. It unions exact Murcko-scaffold groups with broad curated opioid chemical lineages so related opioid chemistry cannot appear in both train and test simply because exact scaffolds differ.

## Fixed baselines

Spectra are represented on nominal integer m/z bins through 700, square-root transformed, and L2-normalized. The models are spectral 1-NN, logistic regression, linear SVM, and a 500-tree random forest with the fixed settings documented in `docs/benchmark_protocol.md`.

| Protocol | Model | Balanced accuracy mean | Balanced accuracy SD | AUROC mean | AUROC SD | Average precision mean |
|---|---|---:|---:|---:|---:|---:|
| Compound-held-out | linear svm | 0.863 | 0.048 | 0.937 | 0.048 | 0.912 |
| Compound-held-out | logistic regression | 0.830 | 0.076 | 0.905 | 0.053 | 0.845 |
| Compound-held-out | random forest | 0.914 | 0.050 | 0.977 | 0.024 | 0.964 |
| Compound-held-out | spectral 1nn | 0.867 | 0.056 | 0.916 | 0.031 | 0.869 |
| Scaffold-held-out | linear svm | 0.734 | 0.151 | 0.751 | 0.270 | 0.700 |
| Scaffold-held-out | logistic regression | 0.709 | 0.147 | 0.696 | 0.322 | 0.651 |
| Scaffold-held-out | random forest | 0.768 | 0.142 | 0.839 | 0.218 | 0.806 |
| Scaffold-held-out | spectral 1nn | 0.722 | 0.155 | 0.707 | 0.331 | 0.694 |
| Analogue-lineage-held-out | linear svm | 0.518 | 0.109 | 0.570 | 0.262 | 0.514 |
| Analogue-lineage-held-out | logistic regression | 0.530 | 0.163 | 0.533 | 0.272 | 0.469 |
| Analogue-lineage-held-out | random forest | 0.581 | 0.116 | 0.685 | 0.206 | 0.597 |
| Analogue-lineage-held-out | spectral 1nn | 0.553 | 0.141 | 0.616 | 0.268 | 0.549 |

## Interpretation

The compound-held-out random forest is strong on this small curated subset (balanced accuracy 0.914; AUROC 0.977), but performance drops under scaffold holdout (0.768; 0.839) and drops further under the stricter analogue-lineage holdout (0.581; 0.685). The analogue result is close enough to chance, with substantial fold-to-fold variation, that novelty performance must not be inferred from the easier compound-held-out result.

## ENFSI external-test feasibility

ENFSI source spectra and row-level derivatives remain local/restricted. The public audit reports aggregate filtering only. Starting from 1,015 ENFSI records that are both metadata-unmatched to SWGDRUG and below the 0.9999 spectral-overlap threshold, conservative ChEBI matching plus three exact-name-and-formula public structure fallbacks resolves 69 structures. Exact InChIKey overlap, exact Murcko-scaffold overlap, and ENFSI internal exact-duplicate filtering leave 23 unique identity/scaffold candidates.

Of those 23 chemically de-leaked candidates, 2 have source-backed opioid labels, 12 have source-backed non-opioid comparator labels, and 9 remain pending. This is not enough positive support for a defensible external classifier metric, so v0.2.0 reports feasibility counts only and deliberately does not report an external AUROC or balanced accuracy.

## Reproducibility and limits

Row-level real-source assignments, predictions, and ENFSI candidate manifests stay local and git-ignored. Public artifacts contain aggregate metrics and hashes only. These results are pilot subset results and do not establish performance on all SWGDRUG compounds, on a broad unseen-opioid distribution, or for operational/confirmatory forensic use.
