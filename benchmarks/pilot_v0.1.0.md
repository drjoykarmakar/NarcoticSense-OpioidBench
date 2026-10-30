# Pilot v0.1.0: leakage-aware EI-MS baselines

This report is an aggregate-only pilot result. Real SWGDRUG split assignments and predictions stay
in ignored local outputs. The supervised subset contains only evidence-curated labels.

## Frozen inputs

- split seed: `20260929`
- folds: 5
- compound-held-out: 184 spectra = 66 opioid + 118 comparator
- scaffold-held-out: 173 spectra = 58 opioid + 115 comparator, 86 exact Murcko scaffold groups
- preprocessing: nominal m/z 0--700, square-root intensity, L2 normalization
- cross-fold exact-spectrum overlap: 0 in both protocols
- cross-fold raw-intensity cosine >= 0.9999: 0 pairs in both protocols
- scaffold overlap in scaffold-held-out CV: 0 in every fold

## Mean 5-fold metrics

| Protocol | Model | Balanced accuracy | AUROC | Average precision | Sensitivity | Specificity |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Compound-held-out | Spectral 1-NN | 0.868 | 0.910 | 0.844 | 0.923 | 0.813 |
| Compound-held-out | Logistic regression | 0.843 | 0.919 | 0.841 | 0.864 | 0.822 |
| Compound-held-out | Linear SVM | 0.880 | 0.964 | 0.939 | 0.878 | 0.881 |
| Compound-held-out | Random forest | 0.913 | 0.981 | 0.974 | 0.834 | 0.991 |
| Scaffold-held-out | Spectral 1-NN | 0.717 | 0.699 | 0.693 | 0.652 | 0.783 |
| Scaffold-held-out | Logistic regression | 0.659 | 0.665 | 0.598 | 0.517 | 0.800 |
| Scaffold-held-out | Linear SVM | 0.729 | 0.752 | 0.725 | 0.614 | 0.843 |
| Scaffold-held-out | Random forest | 0.760 | 0.846 | 0.808 | 0.547 | 0.974 |

Fold-to-fold variance is substantially larger for the scaffold-held-out protocol, consistent with
the small number of positive scaffold groups. These values are preliminary pilot results and must
not be generalized to the 3,801-compound source inventory or to external/unseen-opioid detection.
